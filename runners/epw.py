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
	def generate(self):
		self.getPhonons()
	def getEpw(self):
		pass
	def getPhonons(self):
		mkcd('phonons')
		#cp(dirname(dirname(config.epw.strip()))+'/examples/diamond/pp/C_3.98148.UPF','.')
		#ground-state self-consistent run
		self.preparePh()
		qs=self.getkpoints()
		mkdir('save')
		mkcd('dirs')		
		for i,q in enumerate(qs):
			mkcd(str(i))
			cp('../../scf.in','.')
			passthru(config.mpirun+str(self.m.cores)+config.pw+" < scf.in > scf.out")
			s=read('../../nscf.in').replace('XQ1',str(q[0])).replace('XQ2',str(q[1])).replace('XQ3',str(q[2]))
			write(s,'nscf.in')
			passthru(config.mpirun+str(self.m.cores)+config.pw+" < nscf.in > nscf.out")
			s=read('../../ph.in').replace('XQ1',str(q[0])).replace('XQ2',str(q[1])).replace('XQ3',str(q[2]))
			write(s,'ph.in')
			
			#compute the dynamical matrix, phonon frequencies and change of potential using		
			passthru(config.mpirun+str(self.m.cores)+config.ph+" < ph.in > ph.out ")
			cp('diam.dvscf','../../save/dvscf.%s'%(str(i)))
			cp('diam.dyn','../../save/dyn.%s'%(str(i)))
			cd('..')
	def preparePh(self):
		m=self.m
		pos='  \n'.join(['%s '%(a.symbol)+m.toString(m.atoms.get_scaled_positions()[i]) for i,a in enumerate(m.atoms)])
		cell='\n  '.join([m.toString(m.atoms.cell[i]) for i in range(3)])
		masses=m.getMassFromLabel(m.elements)
		potentials=['C_3.98148.UPF']
		atomspecies='\n'.join(toString(a) for a in zip(m.atoms.get_chemical_symbols(),masses,potentials))
		path=dirname(dirname(config.epw.strip()))+'/examples/diamond/pp/'
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
    celldm(1)       =1.0
    nat             = %d
    ntyp            = %d
    ecutwfc         = 60
    occupations     = 'smearing'
    smearing        = 'mp'
    degauss         = 0.02
    !nbnd            = 4
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
K_POINTS automatic
%s 1 1 1
CELL_PARAMETERS angstrom
%s

"""%(path,len(m.atoms),len(m.elements),atomspecies,pos,toString(m.kpoints),cell)
		write(s.replace('TYPE','scf').replace('PHONON',''),'scf.in')
		write(s.replace('TYPE','phonon').replace('PHONON',"""&phonon
  xqq(1) = XQ1
  xqq(2) = XQ2
  xqq(3) = XQ3
 /
 """),'nscf.in')
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
		write(s,'ph.in')
		q=self.getkpoints()
		np.savetxt('klist.txt',np.array(q))

