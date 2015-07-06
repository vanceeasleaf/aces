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
		cmd=config.phono3py+"--cf3 "
		for file in files:
			dir="dirs/dir_"+file
			cmd+=dir+'/vasprun.xml '
		#generate FORCE_SETS
		passthru(cmd)
		m=self.m
		#Create fc2.hdf and fc3.hdf
		dim=' '.join(map(str,m.supercell))
		passthru(config.phono3py+" --dim='%s'"%(dim))
		
	def generate_supercells(self):
		m=self.m
		#generate supercells
		dim=' '.join(map(str,m.supercell))
		passthru(config.phono3py+"-d --dim='%s'"%(dim))

		
	def generate(self):
		m=self.m
		self.minimizePOSCAR()
		self.generate_supercells()
		files=shell_exec("ls *-*").split('\n')
		
		
		maindir=shell_exec('pwd')
		for file in files:
			print file
			dir="dirs/dir_"+file
			mkdir(dir)
			mv(file,dir+'/POSCAR')
			cd(dir)
			if m.engine=="lammps":			
				self.getVaspRun_lammps()
			elif m.engine=="vasp":
				self.getVaspRun_vasp()
			cd(maindir)
		self.force_constant(files)
		
		#Thermal conductivity calculation
		dim=' '.join(map(str,m.supercell))
		passthru(config.phono3py+'--fc3 --fc2 --dim="'+dim+'" -v --mesh="'+' '.join(map(str,m.kpoints))+'"  --br --thm')
		filename="kappa-m%s.hdf5"%''.join(map(str,m.kpoints))
		passthru(config.pypath+'kaccum --mesh="'+' '.join(map(str,m.kpoints))+'"   '+filename+' |tee kaccum.dat')
		

