#encoding:utf8
from aces.runners import Runner
from ase.io.vasp import write_vasp
from ase import io
from aces.tools import *
import numpy as np
from  aces import config
from aces.runners.phonopy import runner as Runner
from aces.f import read_forces,matrixFormat
import vasp2boltz
from ase.lattice.spacegroup import Spacegroup
class runner(Runner):
	def generate(self):
		cp('minimize/CONTCAR','POSCAR')
		self.prepareBand()
		self.prepareET()
		self.runET()
	def prepareBand(self):
		mkcd('band')
		#run vasp
		cp('../POSCAR','.')
		self.getVaspRun_vasp()
		#cp('/home1/xggong/zhouy/tcscripts/sis/sis1.1/0/secondorder/dirs/dir_POSCAR-002/EIGENVAL','.')
		#cp('/home1/xggong/zhouy/tcscripts/sis/sis1.1/0/secondorder/dirs/dir_POSCAR-002/OUTCAR','.')
		cd('..')
	def prepareET(self):
		atoms = io.read('POSCAR')
		from pyspglib import spglib
		s=spglib.get_spacegroup(atoms)
		from aces.scanf import sscanf
		sgN=sscanf(s,'%s (%d)')[1]
		sg = Spacegroup(sgN)
		atoms.info = {'spacegroup': sg}
		nelect=self.getNelect(outcar='band/OUTCAR')
		mkdir('ET')
		# The remaining part takes care of writing the files required by BoltzTraP.
		bandstructure= vasp2boltz.get_vasp_bandstructure(pathname='band/')
		cd('ET')
		vasp2boltz.write_bandstructure_boltztrap(bandstructure,filename='ET.energy')
		vasp2boltz.write_structure_boltztrap(atoms,filename='ET.struct')
		vasp2boltz.write_intrans_boltztrap(n_electrons = nelect,filename='ET.intrans')
		s=read('ET.energy')
		s=s.replace('HTE','ET')
		write(s,'ET.energy')
		s=read('ET.struct')
		s=s.replace('HTE','ET')
		write(s,'ET.struct')
		cd('..')
	def runET(self):
		cd('ET')
		passthru(config.x_trans+" BoltzTraP")
		cd('..')
	def getNelect(self,outcar='OUTCAR'):
		# Read number of electrons from OUTCAR.
		foundnelect = False
		try:
			for line in open(outcar, 'r'):
				if 'NELECT' in line:
					nelect = float(line.split()[2])
				foundnelect = True
		except:
			pass
		if not foundnelect:
			print 'Number of electrons not found. Please set the number manually in hte.intrans.'
			nelect = 1.0
		return nelect