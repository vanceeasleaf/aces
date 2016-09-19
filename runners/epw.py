#encoding:utf8
from aces.runners import Runner
from ase.io.vasp import write_vasp
from ase import io
from aces.tools import *
import numpy as np
from  aces import config
from aces.runners.phonopy import runner as Runner
from aces.f import read_forces,matrixFormat
from aces.graph import plot,series
class runner(Runner):
	def getNscfKpoints(self):
		c=self.m.ekpoints
		q=[]
		u=c
		for i in range(u[0]):
			for j in range(u[1]):
				for k in range(u[2]):
					b=np.array([float(i)/c[0],float(j)/c[1],float(k)/c[2],1.0/(c[0]*c[1]*c[2])])
					q.append(b)
		return q
	def a(self):
		print config.espresso
	def generate(self):
		self.prepare()
		self.getPhonons()
		self.getEpwIn()
		self.getepw()
	def lifetime(self):
		xx=np.loadtxt('linewidth.elself',skiprows=2)
		E=xx[:,2]
		a=E.argsort()
		h=4.13566743e-15 #eV s	
		h=h*1000*1e12 #meV ps	
		pi2=3.14159265359*2.0
		hbar=h/pi2
		linewidth=xx[:,3] #meV
		lifetime=hbar/linewidth
		plot((E[a],'Electron Energy (meV)'),(lifetime[a],'Electron Lifetime (ps)'),filename='linewidth_elself.png')
	def srunepw1(self):
		m=self.m
		subcores=m.cores
		epw=config.mpirun+str(subcores)+config.epw+" -npool "+str(subcores)
		cmd3=epw+" <epw1.in > epw1.out"
		passthru(cmd3)
	def srunepw(self):
		m=self.m
		subcores=m.cores
		epw=config.mpirun+str(subcores)+config.epw+" -npool "+str(subcores)
		cmd3=epw+" <epw.in > epw.out"
		passthru(cmd3)
	def getepw(self):
		m=self.m
		maindir=pwd()
		cd('epw')
		subcores=m.cores
		pw=config.mpirun+str(subcores)+config.pw+" -npool "+str(subcores)
		epw=config.mpirun+str(subcores)+config.epw+" -npool "+str(subcores)
		cmd1=pw+" <scf.in > scf.out"
		cmd2=pw+" <nscf.in > nscf.out"
		cmd3=epw+" <epw.in > epw.out"
		passthru(cmd1)
		passthru(cmd2)
		passthru(cmd3)
		cd(maindir)
	def collect(self):
		nqpt=len(ls("diam.dyn*"))-1
		for iqpt in np.arange(1,nqpt+1):
			label = str(iqpt)
			shell_exec('cp diam.dyn'+str(iqpt)+' save/diam.dyn_q'+label)
			if (iqpt == 1):
				shell_exec('cp _ph0/diam.dvscf1 save/diam.dvscf_q'+label)
				shell_exec('cp -r _ph0/diam.phsave save/')
			else:
				shell_exec('cp _ph0/diam.q_'+str(iqpt)+'/diam.dvscf1 save/diam.dvscf_q'+label)
				shell_exec('rm _ph0/diam.q_'+str(iqpt)+'/*wfc*' )
	def getPhonons(self):
		m=self.m
		maindir=pwd()
		cd('phonons')
		mkdir('save')
		subcores=m.cores
		pw=config.mpirun+str(subcores)+config.pw+" -npool "+str(subcores)
		ph=config.mpirun+str(subcores)+config.ph+" -npool "+str(subcores)
		cmd1=pw+" < scf.in  > scf.out"
		cmd2=ph+" < ph.in   > ph.out "
		passthru(cmd1)	
		passthru(cmd2)
		self.collect()
		cd(maindir)
	def getPotName(self):
		maindir=pwd()
		cd(config.qepot)
		files=ls("*")
		cd(maindir)
		import re  
		q=[re.split("[_\.]",a)[0].capitalize() for a in files]
		u={}
		for i,v in enumerate(q):
			u[v]=files[i]
		return u
	def get_amass(self):
		m=self.m
		masses=m.getMassFromLabel(m.elements)
		amass='\n'.join(['amass(%d)=%f'%(i+1,s) for i,s in enumerate(masses)])
		return amass
	def get_dos(self):
		m=self.m
		amass=self.get_amass()
		nks="nk1=15,nk2= 15,nk3= 15"
		write("""&input
    asr='simple',  
    dos=.true. 
    %s
    flfrc='diam.fc', 
    fldos='diam.dos', 
    %s
 /"""%(amass,nks),"phdos.in")
		shell_exec(config.espresso+"bin/matdyn.x < phdos.in >phdos.out")
		self.drawDos()
	def readdos(self):
		xx=np.loadtxt('diam.dos')
		freq=xx[:,0]
		dos=xx[:,1]
		return freq,dos
	def getbanddos(self):
		freq,dos=self.readdos()
		freq/=33.367
		from aces.bandplot import plotbanddos
		plotbanddos(freq,dos,labels=' '.join(self.m.bandpath))
	def drawDos(self):
		freq,dos=self.readdos()
		freq/=33.367
		plot((freq,'Frequency (THz)'),(dos,'Density of States'),filename='total_dos.png')
	def get_band(self):
		
		amass=self.get_amass()
		m=self.m
		bp=m.bandpoints
		rcell=m.atoms.get_reciprocal_cell().T
		kp =[rcell.dot(bp[x]) for x in m.bandpath]
		bpath=[m.toString(x) for x in kp]
		bs='\t101\n'.join(bpath)
		write("""&input
    asr='simple',  
    %s
    flfrc='diam.fc', 
    flfrq='diam.freq', 
    q_in_band_form=.true.,
 /
  %d
  %s 1"""%(amass,len(m.bandpath),bs),"band.in")
		shell_exec(config.espresso+"bin/matdyn.x < band.in >band.out")

		self.get_bandplot()
		
		a=np.loadtxt("band.plot")
		c=np.unique(a[:,0])
		c=c[0:-1]
		nqpoint=len(c)
		s="nqpoint: %d\n"%(nqpoint)
		s+="""npath: %d\n"""%(len(self.m.bandpath)-1)
		s+="phonon: \n"
		for x in c:
			filter=(a[:,0]==x)
			ys=a[filter,1]
			s+="- q-position: [    0,    0.0000000,    0.0000000 ]\n"
			s+="  distance:    %f\n"%x
			s+="  band:\n"
			for i,y in enumerate(ys):
				s+="""  - # %d
    frequency:   %f\n"""%(i+1,y/33.367)
			s+="\n"
		write(s,"band.yaml")
		
		
	def get_bandplot(self):
		freq,dos=self.readdos()
		write("""diam.freq
0 %f
band.plot
band.ps
0.0
50.0 0.0"""%(freq.max()),"plotband.in")
		shell_exec(config.espresso+"bin/plotband.x < plotband.in >plotband.out")
	def post(self):
		self.get_force_constants()
		self.get_dos()
		self.get_band()		
		self.getbanddos()
	def get_force_constants(self):
		write("""&input
   fildyn='diam.dyn', zasr='simple', flfrc='diam.fc'
 /""","q2r.in")
		shell_exec(config.espresso+"bin/q2r.x < q2r.in >q2r.out")
	def getKlist(self,dir="./"):
		q=np.loadtxt(dir+"diam.dyn0",skiprows=2)
		assert len(q)>0
		return q
		#another way ,which maybe not acurrate and cause some bugs
		file=open(dir+"ph.out")
		q=[]
		from aces.scanf import sscanf
		for line in file:
			if "       N         xq(1)         xq(2)         xq(3)" in line:
				while True:
					line=file.next().strip()
					if(line==""):
						break
					k=sscanf(line,"%d   %f   %f   %f")[1:]
					q.append(k)
				break
		assert len(q)>0
		return q
	def getEpwIn(self):
		m=self.m
		nks="nk1         = %d\n  nk2         = %d\n  nk3         = %d"%tuple(m.ekpoints)
		nkfs="nkf1         = %d\n  nkf2         = %d\n  nkf3         = %d"%tuple(m.nkf)
		nqs="nq1         = %d\n  nq2         = %d\n  nq3         = %d"%tuple(m.kpoints)
		nqfs="nqf1         = %d\n  nqf2         = %d\n  nqf3         = %d"%tuple(m.nqf)
		amass=self.get_amass()
		q=self.getKlist(dir="phonons/")
		qs='\n'.join([m.toString(a)+" %f"%(1.0/len(q)) for a in q])
		s="""--
&inputepw
  prefix      = 'diam'
  %s
  outdir      = './'

  elph        = .true.
  kmaps       = .false.
  epbwrite    = .true.
  epbread     = .false.

  epwwrite    = .true.
  epwread     = .false.

  nbndskip    =  0

  wannierize  = .true.
  num_iter    = 300
  iprint      = 2
  dis_win_max = 12
  dis_froz_max= 7
  proj(1)     = 'random'   

  iverbosity  = 0

  elinterp    = .true.
  phinterp    = .true.

  tshuffle2   = .true.
  tphases     = .false.

  elecselfen  = .true.
  phonselfen  = .true.
  a2f         = .true.

  parallel_k  = .true.
  parallel_q  = .false.

  fsthick     = 2.80284905 ! eV
  eptemp      = 300 ! K
  degaussw    = 0.1 ! eV

  dvscf_dir   = '../phonons/save'
  filukk      = './diam.ukk'
  
  %s
  %s
  %s
  %s

 /
  %d cartesian
%s

"""%(amass,nqs,nks,nqfs,nkfs,len(q),qs)

		write(s,"epw/epw.in")

	def prepare(self):
		"""为了保证q+k也在k的格点上，必须要求qs*cell*[nk1,nk2,nk3]为整数，即nk=N*nq
		否则可能出现Error in routine createkmap (1):
		"""
		m=self.m
		u=np.array(m.ekpoints)/np.array(m.kpoints)
		assert np.allclose(np.floor(u+.5),u,rtol=0.01)
		m=self.m
		mkdir("phonons")
		mkdir("epw")
		pos='  \n'.join(['%s '%(a.symbol)+m.toString(m.atoms.get_scaled_positions()[i]) for i,a in enumerate(m.atoms)])
		mmm=1.0#np.abs(m.atoms.cell).max()*2.0
		cell='\n  '.join([m.toString(m.atoms.cell[i]/mmm) for i in range(3)])
		masses=m.getMassFromLabel(m.elements)
		pots=self.getPotName()
		potentials=[pots[a] for a in m.elements]
		atomspecies='\n'.join(toString(a) for a in zip(m.elements,masses,potentials))
		path=config.qepot
		skpoints="K_POINTS automatic\n%s 1 1 1"%toString(m.ekpoints)
		s=""" &control
    calculation     = 'TYPE'
    prefix          = 'diam'
    restart_mode    = 'from_scratch'
    wf_collect      = .false.
    pseudo_dir      = '%s'
    outdir          = './'
    tprnfor         = .true.
    tstress         = .true.
    etot_conv_thr = 1.0d-5
    forc_conv_thr = 1.0d-4
 /
 &system
    ibrav           = 0
    celldm(1)       =%f
    nat             = %d
    ntyp            = %d
    ecutwfc         = 40
    occupations     = 'smearing'
    smearing        = 'mp'
    degauss         = 0.02
 /
 &electrons
    diagonalization = 'david'
    mixing_beta     = 0.7
    conv_thr        = 1.0d-8
 /
PHONON
ATOMIC_SPECIES
  %s
ATOMIC_POSITIONS crystal
%s
%s
CELL_PARAMETERS alat
%s

"""%(path,mmm/0.5291772083,len(m.atoms),len(m.elements),atomspecies,pos,skpoints,cell)
		write(s.replace('TYPE','scf').replace('PHONON',''),'phonons/scf.in')
		cp('phonons/scf.in',"epw/")

		nscfk=self.getNscfKpoints()
		aa='\n'.join([m.toString(a) for a in nscfk])
		newskpoints="K_POINTS crystal\n%d\n%s"%(len(nscfk),aa)
		write(
			s.replace('TYPE','nscf')
			.replace("restart_mode    = 'from_scratch'","")
			.replace("tprnfor         = .true.","")
			.replace("tstress         = .true.","")
			.replace(skpoints,newskpoints)
			,
		"epw/nscf.in")
		amass=self.get_amass()
		nqs="nq1=%d,\n  nq2= %d,\n  nq3= %d,"%tuple(m.kpoints)
		s="""--
&inputph
  prefix   = 'diam'
  fildyn   = 'diam.dyn'
  fildvscf = 'dvscf'
  ldisp    = .true
  %s
  tr2_ph   =  1.0d-10
 /
"""%nqs
		write(s,'phonons/ph.in')

