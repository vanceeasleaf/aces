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
		self.getPhonons()
		self.getepw()
		self.getbol()
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
		plot((E[a],'Electron Energy (eV)'),(lifetime[a],'Electron Lifetime (ps)'),filename='linewidth_elself.png')
		hist, bin_edges = np.histogram(E, density=True,bins=200)
		plot((bin_edges[:-1],'Electron Energy (eV)'),(hist,'Density of States'),filename='fine_edos.png')
	def runepw1(self):
		m=self.m
		subcores=m.cores
		epw=config.mpirun+str(subcores)+config.epw+" -npool "+str(subcores)
		cmd3=epw+" <epw1.in > epw1.out"
		passthru(cmd3)
	def rerunepw(self):
		m=self.m
		subcores=m.cores
		if not exists('epw.in'):
			self.getepwin()
		epw=config.mpirun+str(subcores)+config.epw+" -npool "+str(subcores)
		cmd3=epw+" <epw.in > epw.out"
		passthru(cmd3)
	def runepw(self):
		self.getepwin()
		self.rerunepw()
	def postbol(self):
		seebeck=np.loadtxt("diam_seebeck.dat",skiprows=2)
		ef=self.getfermi()
		mu=seebeck[:,0]-ef
		filter=(mu<2) * (mu>-2)
		mu=mu[filter]
		seebeck=seebeck[filter]
		s1=seebeck[:,2]
		s2=seebeck[:,6]
		s3=seebeck[:,10]
		from aces.graph import fig,pl
		with fig("Seebeck.png",legend=True):
			pl.plot(mu,s1,label="Sxx")
			pl.plot(mu,s2,label="Syy")
			pl.plot(mu,s3,label="Szz")
			pl.xlabel("$\mu$-Ef (eV)")
			pl.ylabel("Seebeck Coefficient (V/K)")
		sigma=np.loadtxt("diam_elcond.dat")
		sigma=sigma[filter]
		with fig("Sigma.png",legend=True):
			pl.plot(mu,sigma[:,2],label="$\sigma$ xx")
			pl.plot(mu,sigma[:,4],label="$\sigma$ yy")
			pl.plot(mu,sigma[:,7],label="$\sigma$ zz")
			pl.xlabel("$\mu$-Ef (eV)")
			pl.ylabel("Electrical conductivity (1/$\Omega$/m)")
		kappa=np.loadtxt("diam_kappa.dat")
		kappa=kappa[filter]
		with fig("kappa.png",legend=True):
			pl.plot(mu,kappa[:,2],label="$\kappa$ xx")
			pl.plot(mu,kappa[:,4],label="$\kappa$ yy")
			pl.plot(mu,kappa[:,7],label="$\kappa$ zz")
			pl.xlabel("$\mu$-Ef (eV)")
			pl.ylabel("Thermal Conductivity (W/m/K)")
		with fig("powerfactor.png",legend=True):
			pl.plot(mu,sigma[:,2]*s1*s1,label="Pxx")
			pl.plot(mu,sigma[:,4]*s2*s2,label="Pyy")
			pl.plot(mu,sigma[:,7]*s3*s3,label="Pzz")
			pl.xlabel("$\mu$-Ef (eV)")
			pl.ylabel("Power Factor ")
		dos=np.loadtxt("diam_boltzdos.dat")
		with fig("boltzdos.png"):
			pl.plot(dos[:,0],dos[:,1])
			pl.xlabel("E-Ef (eV)")
			pl.ylabel("Electron Density of States (arbitrary unit)")
	def getwin(self):
		m=self.m
		pos='  \n'.join(['%s '%(a.symbol)+m.toString(m.atoms.get_scaled_positions()[i]) for i,a in enumerate(m.atoms)])
		cell='\n  '.join([m.toString(x) for x in m.atoms.cell])
		nks=m.toString(m.ekpoints)
		nkfs=m.toString(m.nkf)
		nscfk=self.getNscfKpoints()
		aa='\n'.join([m.toString(a[:-1]) for a in nscfk])
		if m.ekpoints[2]==1:
			boltz_2d_dir='z'
		else:
			boltz_2d_dir='no'
		boltz_2d_dir='z'
		s="""!!! -- Begin of BoltzWann input -- !!!
boltzwann                    = true
boltz_calc_also_dos          = true
boltz_dos_energy_step        = 0.01
boltz_2d_dir                 = %s
smr_type                     = gauss
boltz_dos_adpt_smr           = false
boltz_dos_smr_fixed_en_width = 0.03
kmesh                        = %s
boltz_mu_min                 = -4.
boltz_mu_max                 = 4.
boltz_mu_step                = 0.01
boltz_temp_min               = 300.
boltz_temp_max               = 300.
boltz_temp_step              = 50
boltz_relax_time             = 10.
!!! --- End of BoltzWann input --- !!!     
num_wann          = 14
search_shells=50
iprint   2
dis_win_max       = 17.d0
dis_froz_max      = 6.4d0
dis_num_iter      = 120
dis_mix_ratio     = 1.d0
num_iter     300
num_print_cycles  = 50
begin unit_cell_cart
%s
end unit_cell_cart
begin atoms_frac
%s
End atoms_frac
begin projections     
random
end projections
mp_grid =     %s   
begin kpoints
%s
end kpoints
		"""%(boltz_2d_dir,nkfs,cell,pos,nks,aa)
		#pw2wannier90 make some operation to cell and finally judge rcell == acell (<1e-6) and it's too strict 
		#so we have to change constant eps6to 1e-5
		#im espresso/Modules/constant.f90
        #search_shells default=12 but I change it 50 to avoid expected
        #termination of wannier90
		write(s,"diam.win")
		s="""&inputpp 
   outdir = './'
   prefix = 'diam'
   seedname='diam'
   spin_component = 'none'
   write_mmn = .true.
   write_amn = .true.
   write_unk = .false.
/"""
		write(s,"diam.pw2wan")
	def runbol(self):
		self.getwin()
		p=config.espresso+"bin/postw90.x";
		if not exists(p):
			print "p not exists"
			return
		postw90=self.getx("postw90",pool=False)
		passthru(postw90+" diam")
	def runw(self):
		self.getwin()
		wannier=self.getx("wannier90",pa=False)
		passthru(wannier+" -pp diam")
		pw2w=self.getx("pw2wannier90",pool=False)
		passthru(pw2w +"<diam.pw2wan >pw2wan.out")
		passthru(wannier+"  diam")
		self.runbol()
	def getbol(self):
		m=self.m
		maindir=pwd()
		mkcd('bol')
		self.runscf()
		self.runnscf()
		self.runw()
		cd(maindir)
	def getepw(self):
		m=self.m
		maindir=pwd()
		mkcd('epw')
		"""为了保证q+k也在k的格点上，必须要求qs*cell*[nk1,nk2,nk3]为整数，即nk=N*nq
		否则可能出现Error in routine createkmap (1):
		"""
		u=np.array(m.ekpoints)/np.array(m.kpoints)
		assert np.allclose(np.floor(u+.5),u,rtol=0.01)
		self.runscf()
		self.runnscf()
		self.runepw()
		cd(maindir)
	def runscf(self):
		if not exists('scf.in'):
			self.getscfin()
		pw=self.getx("pw")
		passthru(pw+" <scf.in > scf.out")
	def runnscf(self):
		if not exists('nscf.in'):
			self.getnscfin()
		pw=self.getx("pw")
		passthru(pw+" <nscf.in > nscf.out")
	def collect(self):
		mkdir('save')
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
	def getphin(self):
		m=self.m
		nqs="nq1=%d,\n  nq2= %d,\n  nq3= %d,"%tuple(m.kpoints)
		s="""--
&inputph
  prefix   = 'diam'
  fildyn   = 'diam.dyn'
  fildvscf = 'dvscf'
  ldisp    = .true
  %s
  tr2_ph   =  1.0d-12
 /
"""%nqs
		write(s,'ph.in')
	def runph(self):
		m=self.m
		subcores=m.cores
		if not exists('ph.in'):
			self.getphin()
		ph=self.getx("ph")
		passthru(ph+" < ph.in   > ph.out ")
	def getPhonons(self):
		m=self.m
		maindir=pwd()
		mkcd('phonons')
		self.runscf()
		self.runph()
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
		nks="nk1=%d,nk2= %d,nk3= %d"%tuple(m.kpoints)
		write("""&input
    asr='simple',  
    dos=.true. 
    %s
    flfrc='diam.fc', 
    fldos='diam.dos', 
    %s
 /"""%(amass,nks),"phdos.in")
		matdyn=self.getx("matdyn")
		shell_exec(matdyn+"< phdos.in >phdos.out")
		self.drawDos()
	def readdos(self):
		xx=np.loadtxt('diam.dos')
		freq=xx[:,0]
		dos=xx[:,1]
		return freq,dos
	def getfermi(self):
		a=shell_exec("grep Fermi nscf.out").strip()
		from aces.scanf import sscanf
		fermi=sscanf(a,"the Fermi energy is    %f ev")[0]
		return fermi
	def getedos(self):
		m=self.m
		fermi=self.getfermi()
		write("""&DOS
outdir='./'
prefix='diam'
fildos='diam.edos',
Emin=%f, Emax=%f, DeltaE=0.01,degauss=0.005,ngauss=1
/
"""%(-4.0+fermi,4.0+fermi),"edos.in")
		dos=self.getx("dos")
		shell_exec(dos +" < edos.in >edos.out")
		self.drawedos()
	def getpedos(self):
		m=self.m
		write("""&projwfc
outdir='./'
prefix='diam'
filpdos='diam.pedos',
Emin=-10.0, Emax=10.0, DeltaE=0.1,degauss=0.0025,ngauss=-99
/
""","pedos.in")
		projwfc=self.getx('projwfc')
		shell_exec(projwfc+" < pedos.in >pedos.out")
	def getx(self,s='pw',pa=True,pool=True):
		subcores=self.m.cores
		if not pa:
			return config.espresso+"bin/"+s+".x "
		s= config.mpirun+str(subcores)+" "+config.espresso+"bin/"+s+".x "
		if pool:
			s+="-npool "+str(subcores)
		return s
	def drawpedos(self):
		s=np.loadtxt("diam.edos",skiprows=1)
		freq,dos=s[:,0],s[:,1]
		plot((freq,'Electron Energy (eV)'),(dos,'Density of States'),filename='edos.png')

	def getbanddos(self):
		freq,dos=self.readdos()
		freq/=33.367
		from aces.bandplot import plotbanddos
		plotbanddos(freq,dos,labels=' '.join(self.m.bandpath))
	def drawedos(self):
		fermi=self.getfermi()
		s=np.loadtxt("diam.edos",skiprows=1)
		freq,dos=s[:,0],s[:,1]
		freq-=fermi
		plot((freq,'Electron Energy (eV)'),(dos,'Density of States'),filename='edos.png')

	def drawDos(self):
		freq,dos=self.readdos()
		freq/=33.367
		plot((freq,'Frequency (THz)'),(dos,'Density of States'),filename='total_dos.png')
	def get_band(self):
		
		amass=self.get_amass()
		m=self.m
		bp=m.bandpoints
		bpath=[m.toString(bp[x]) for x in m.bandpath]
		bs='\t100\n'.join(bpath)
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
		write("""diam.freq
0 600
band.plot
band.ps
0.0
50.0 0.0""","plotband.in")
		shell_exec(config.espresso+"bin/plotband.x < plotband.in >plotband.out")
	def post(self):
		self.get_force_constants()
		self.get_band()
		self.get_dos()
		self.getbanddos()
	def get_force_constants(self):
		write("""&input
   fildyn='diam.dyn', zasr='simple', flfrc='diam.fc'
 /""","q2r.in")
		shell_exec(config.espresso+"bin/q2r.x < q2r.in >q2r.out")
	def getKlist(self,dir="./"):
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
	def getepwin(self):
		m=self.m
		nks="nk1         = %d\n  nk2         = %d\n  nk3         = %d"%tuple(m.ekpoints)
		nkfs="nkf1         = %d\n  nkf2         = %d\n  nkf3         = %d"%tuple(m.nkf)
		nqs="nq1         = %d\n  nq2         = %d\n  nq3         = %d"%tuple(m.kpoints)
		nqfs="nqf1         = %d\n  nqf2         = %d\n  nqf3         = %d"%tuple(m.nqf)
		amass=self.get_amass()
		q=self.getKlist(dir="../phonons/")
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
  wdata(1)    = 'search_shells=50'

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
		write(s,"epw.in")
	def getscf_template(self):
		m=self.m
		pos='  \n'.join(['%s '%(a.symbol)+m.toString(m.atoms.get_scaled_positions()[i]) for i,a in enumerate(m.atoms)])
		mmm=1.0#np.abs(m.atoms.cell).max()*2.0
		cell='\n  '.join([m.toString(m.atoms.cell[i]/mmm) for i in range(3)])
		masses=m.getMassFromLabel(m.elements)
		pots=self.getPotName()
		potentials=[pots[a] for a in m.elements]
		atomspecies='\n'.join(toString(a) for a in zip(m.elements,masses,potentials))
		path=config.qepot
		
		s=""" &control
    calculation     = 'TYPE'
    prefix          = 'diam'   
    wf_collect      = .true.
    pseudo_dir      = '%s'
    outdir          = './'   
    CONTROL
    etot_conv_thr = 1.0d-5
    forc_conv_thr = 1.0d-4
 /
 &system
    ibrav           = 0
    celldm(1)       =%f
    nat             = %d
    ntyp            = %d
    ecutwfc         = 60
    occupations     = 'smearing'
    smearing        = 'mp'
    degauss         = 0.02
 /
 &electrons
    diagonalization = 'david'
    mixing_beta     = 0.7
    conv_thr        = 1.0d-10
 /
ATOMIC_SPECIES
  %s
ATOMIC_POSITIONS crystal
%s
KPOINTS
CELL_PARAMETERS alat
%s

"""%(path,mmm/0.5291772083,len(m.atoms),len(m.elements),atomspecies,pos,cell)
		return s
	def getscfin(self):
		m=self.m
		tmpl=self.getscf_template()
		skpoints="K_POINTS automatic\n%s 1 1 1"%toString(m.ekpoints)
		a="""tprnfor         = .true.
    tstress         = .true.
    restart_mode    = 'from_scratch'
    """
		s=tmpl.replace('TYPE','scf').replace('CONTROL',a).replace('KPOINTS',skpoints)
		write(s,'scf.in')
	def getnscfin(self):
		m=self.m
		tmpl=self.getscf_template()
		nscfk=self.getNscfKpoints()
		aa='\n'.join([m.toString(a) for a in nscfk])
		newskpoints="K_POINTS crystal\n%d\n%s"%(len(nscfk),aa)
		s=tmpl.replace('TYPE','nscf').replace('CONTROL','').replace('KPOINTS',newskpoints)
		write(s,"nscf.in")

