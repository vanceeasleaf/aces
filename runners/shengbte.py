#encoding:utf8
from aces.tools import *
import aces.config as config
from ase.io import read
from ase.io.vasp import write_vasp
from aces.binary import pr
from aces.runners import Runner
from aces.UnitCell.unitcell import UnitCell
from aces.graph import plot,series
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pl
import numpy as np
from aces.runners.phonopy import runner as Runner
class runner(Runner):
			
	def force_constant(self,files):
		cmd=config.phonopy+"-f "
		for file in files:
			dir="dirs/dir_"+file
			cmd+=dir+'/vasprun.xml '
		#generate FORCE_SETS
		passthru(cmd)
		m=self.m
		#Create FORCE_CONSTANTS
		passthru(config.phonopy+"--writefc --dim='%s'"%(m.dim))
	
	def force_constant3(self,files):
		m=self.m
		cmd='find dirs/dir_3RD.* -name vasprun.xml |sort -n|'+config.thirdorder+" reap %s 0.54 "%m.dim
		passthru(cmd)

	def generate_supercells3(self):
		m=self.m
		#generate supercells
		passthru(config.thirdorder+"sow %s 0.54"%(m.dim))
	
	def getControl(self):
		m=self.m
		f=open('CONTROL','w')
		atoms=read('../POSCAR')#m.atoms
		allocations="""&allocations
	nelements=%d
	natoms=%d
	ngrid(:)=%s
&end
"""%(len(m.elements),len(atoms),' '.join(map(str,m.kpoints)))
		cell=atoms.cell
		types=m.toString([m.elements.index(x)+1 for x in atoms.get_chemical_symbols()])
		pos=""
		for i,atom in enumerate(atoms):
			pos+="	positions(:,%d)=%s\n"%(i+1,m.toString(atoms.get_scaled_positions()[i]))
		crystal="""&crystal
	lfactor=0.1,
	lattvec(:,1)=%s
	lattvec(:,2)=%s
	lattvec(:,3)=%s
	elements=%s
	types=%s
%s
	scell(:)=%s
&end
"""%(' '.join(map(str,cell[0])),
			' '.join(map(str,cell[1])),
			' '.join(map(str,cell[2])),
			' '.join(map(lambda x: '"'+x+'"',m.elements)),
			types,
			pos,
			m.dim)
		parameters="""&parameters
	T=%f
	scalebroad=1.0
&end
"""%(m.T)
		flags="""
&flags
	nonanalytic=.TRUE.
	nanowires=.FALSE.
&end  
		"""
		f.write(allocations)
		f.write(crystal)
		f.write(parameters)
		f.write(flags)
		f.close()

	def generate(self):
		m=self.m
		self.minimizePOSCAR()
		#cp('minimize/POSCAR','.')
		mkdir('secondorder')
		cd('secondorder')
		cp('../POSCAR','.')
		self.generate_supercells()
		files=shell_exec("ls *-*").split('\n')
		assert len(files)>0
		self.getvasprun(files)
		self.force_constant(files)
		cd('..')
		mkdir('thirdorder')
		cd('thirdorder')
		cp('../POSCAR','.')
		self.generate_supercells3()
		files=shell_exec("ls 3RD.*.*|sort -n").split('\n')
		assert len(files)>0
		self.getvasprun(files)
		self.force_constant3(files)
		cd('..')
		mkdir('SHENG')
		cd('SHENG')
		cp('../secondorder/FORCE_CONSTANTS','FORCE_CONSTANTS_2ND')
		cp('../thirdorder/FORCE_CONSTANTS_3RD','.')
		self.getControl()
		#Thermal conductivity calculation
		print "START SHENGBTE..."
		passthru(config.mpirun+" %s "%m.cores+config.shengbte)
		

