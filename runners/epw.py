#encoding:utf8
from aces.runners import Runner
from ase.io.vasp import write_vasp
from ase import io
from aces.tools import *
import numpy as np
from  aces import config
from aces.runners.phonopy import runner as Runner
from aces.f import read_forces,matrixFormat
class runner(Runner):
	def getkpoints(self):
		c=self.m.kpoints
		q=[]
		u=[int(x/2)*2+1 for x in c]
		for i in range(u[0]):
			for j in range(u[1]):
				for k in range(u[2]):
					b=np.array([float(i-c[0]/2)/c[0],float(j-c[1]/2)/c[1],float(k-c[2]/2)/c[2]])
					q.append(b)
		return q
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
	def generate(self):
		self.prepare()
		self.getPhonons()
		self.getepw()
	def getepw(self):
		m=self.m
		maindir=pwd()
		cd('epw')
		mkdir('out')
		mkcd('tmp')
		subcores=4
		pw=config.mpirun+str(subcores)+config.pw+" -npool "+str(subcores)
		epw=config.mpirun+str(subcores)+config.epw+" -npool "+str(subcores)
		cmd1=pw+" <../inp/scf.in > ../out/scf.out"
		cmd2=pw+" <../inp/nscf.in > ../out/nscf.out"
		cmd3=epw+" <../inp/epw.in > ../out/epw.out"
		from aces.jobManager import jobManager
		self.jm=jobManager()
		if 'jm' in self.__dict__:

			from aces.jobManager import pbs
			path=pwd()
			if m.queue=="q3.4":
				pb=pbs(queue=m.queue,nodes=subcores,procs=1,disp=m.pbsname,path=path,content=cmd1+"\n"+cmd2+"\n"+cmd3)
			else:
				pb=pbs(queue=m.queue,nodes=1,procs=subcores,disp=m.pbsname,path=path,content=cmd1+"\n"+cmd2+"\n"+cmd3)
			self.jm.reg(pb)	
		self.jm.check()
		#passthru(cmd1)
		#passthru(cmd2)
		#passthru(cmd3)
	def getPhonons(self):
		m=self.m
		maindir=pwd()
		cd('phonons')
		qs=self.getkpoints()
		mkdir('save')
		mkcd('tmp')		
		subcores=4
		pw=config.mpirun+str(subcores)+config.pw+" -npool "+str(subcores)
		ph=config.mpirun+str(subcores)+config.ph+" -npool "+str(subcores)
		if subcores < 10:
			dvscf="diam.dvscf1"
		else:
			dvscf="diam.dvscf01"
		from aces.jobManager import jobManager
		self.jm=jobManager()
		for i,q in enumerate(qs):
			mkcd(str(i))
			cp('../../inp/scf.in','.')					
			s=read('../../inp/nscf.in').replace('XQ1',str(q[0])).replace('XQ2',str(q[1])).replace('XQ3',str(q[2]))
			write(s,'nscf.in')			
			s=read('../../inp/ph.in').replace('XQ1',str(q[0])).replace('XQ2',str(q[1])).replace('XQ3',str(q[2]))
			write(s,'ph.in')
			mkcd("tmp")	
			cmd1=pw+" < ../scf.in  > ../scf.out"
			cmd2=pw+" < ../nscf.in > ../nscf.out"
			cmd3=ph+" < ../ph.in   > ../ph.out "
			cmd4="cp "+dvscf+" ../../../save/dvscf.%s"%(str(i))+";cp diam.dyn ../../../save/dyn.%s"%(str(i))
			if 'jm' in self.__dict__:
				from aces.jobManager import pbs
				path=pwd()
				if m.queue=="q3.4":
					pb=pbs(queue=m.queue,nodes=subcores,procs=1,disp=m.pbsname,path=path,content=cmd1+"\n"+cmd2+"\n"+cmd3+"\n"+cmd4)
				else:
					pb=pbs(queue=m.queue,nodes=1,procs=subcores,disp=m.pbsname,path=path,content=cmd1+"\n"+cmd2+"\n"+cmd3+"\n"+cmd4)
				self.jm.reg(pb)			
			else:	
				passthru(cmd1)	
				passthru(cmd2)
				passthru(cmd3)
				passthru(cmd4)
			cd('../..')
		self.jm.check()
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
	def prepare(self):
		m=self.m
		mkdir("phonons/inp")
		mkdir("epw/inp")
		pos='  \n'.join(['%s '%(a.symbol)+m.toString(m.atoms.get_scaled_positions()[i]) for i,a in enumerate(m.atoms)])
		mmm=np.abs(m.atoms.cell).max()*2.0
		cell='\n  '.join([m.toString(m.atoms.cell[i]/mmm) for i in range(3)])
		masses=m.getMassFromLabel(m.elements)
		pots=self.getPotName()
		potentials=[post[a] for a in m.elements]
		atomspecies='\n'.join(toString(a) for a in zip(m.atoms.get_chemical_symbols(),masses,potentials))
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
    nbnd            = 4
 /
 &electrons
    diagonalization = 'david'
    mixing_beta     = 0.7
    conv_thr        = 1.0d-10
 /
PHONON
ATOMIC_SPECIES
  %s
ATOMIC_POSITIONS alat
%s
%s
CELL_PARAMETERS angstrom
%s

"""%(path,mmm/0.5291772083,len(m.atoms),len(m.elements),atomspecies,pos,skpoints,cell)
		write(s.replace('TYPE','scf').replace('PHONON',''),'phonons/inp/scf.in')
		cp('phonons/inp/scf.in',"epw/inp/")
		write(s.replace('TYPE','phonon').replace('PHONON',"""&phonon
  xqq(1) = XQ1
  xqq(2) = XQ2
  xqq(3) = XQ3
 /
 """),'phonons/inp/nscf.in')
		nscfk=self.getNscfKpoints()
		aa='\n'.join([m.toString(a) for a in nscfk])
		newskpoints="K_POINTS crystal\n%d\n%s"%(len(nscfk),aa)
		write(
			s.replace('TYPE','nscf').replace("restart_mode    = 'from_scratch'","").replace("tprnfor         = .true.","").replace("tstress         = .true.","")
			.replace(skpoints,newskpoints)
			,"epw/inp/nscf.in")
		amass='\n'.join(['amass(%d)=%f'%(i+1,s) for i,s in enumerate(masses)])
		s="""--
&inputph
  tr2_ph   =  1.0d-12
  prefix   = 'diam'
  %s
  outdir   = './'
  fildyn   = 'diam.dyn'
  fildvscf = 'dvscf'
 /
XQ1 XQ2 XQ3
"""%amass
		write(s,'phonons/inp/ph.in')
		q=self.getkpoints()
		np.savetxt('phonons/klist.txt',np.array(q))
		nks="nk1         = %d\n  nk2         = %d\n  nk3         = %d"%tuple(m.ekpoints)
		nqs="nq1         = %d\n  nq2         = %d\n  nq3         = %d"%tuple(m.kpoints)
		nkfs="nkf1         = %d\n  nkf2         = %d\n  nkf3         = %d"%tuple(m.ekpoints)
		nqfs="nqf1         = %d\n  nqf2         = %d\n  nqf3         = %d"%tuple(m.kpoints)
		qs=' .5\n'.join([m.toString(a) for a in q])
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

  epf_mem     = .true.

  nbndsub     =  4
  nbndskip    =  0

  wannierize  = .true.
  num_iter    = 300
  iprint      = 2
  dis_win_max = 12
  dis_froz_max= 7
  proj(1)     = 'Si:sp3'   

  iverbosity  = 0

  elinterp    = .true.
  phinterp    = .true.

  tshuffle2   = .true.
  tphases     = .false.

  elecselfen  = .false.
  phonselfen  = .true.
  a2f         = .true.

  parallel_k  = .true.
  parallel_q  = .false.

  fsthick     = 6.80284905 ! eV
  eptemp      = 0.075 ! K
  degaussw    = 1.0204273575 ! eV

  dvscf_dir   = '../../phonons/save'
  filukk      = './SiC.ukk'
  
  %s
  %s
  %s
  %s

! The list of qpoints below must be identical
! to ../phonons/qlist.dat
 /
  %d cartesian
%s
"""%(amass,nqs,nks,nqfs,nkfs,len(q),qs)
		write(s,"epw/inp/epw.in")
