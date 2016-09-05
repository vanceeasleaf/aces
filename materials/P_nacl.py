from aces.material import material
from aces.modify import get_unique_atoms
from ase import Atoms,Atom
from math import pi,sqrt
from ase.dft.kpoints import ibz_points
from aces import config
from ase.lattice import bulk
import numpy as np
from aces.tools import *
class structure(material):
	def set_parameters(self):
		self.modulation=False
		pass#['Gamma','Y','T','X','Gamma']

	def setup(self):
		self.forceThick=False
		self.elements=['P']
		self.bandpoints=ibz_points['fcc']
		self.bandpath=['Gamma','K','X','Gamma','L']
		self.premitive/=np.array([self.latx,self.laty,self.latz])

	def lmp_structure(self):
		pos=np.array([
			[  0.0000000000000000, 0.0000000000000000,  0.0000000000000000],
			[  0.0000000000000000, 0.5000000000000000,  0.5000000000000000],
			[  0.5000000000000000, 0.0000000000000000,  0.5000000000000000],
			[  0.5000000000000000, 0.5000000000000000,  0.0000000000000000],
			[  0.5000000000000000, 0.5000000000000000,  0.5000000000000000],
			[  0.5000000000000000, 0.0000000000000000,  0.0000000000000000],
			[  0.0000000000000000, 0.5000000000000000,  0.0000000000000000],
			[  0.0000000000000000, 0.0000000000000000,  0.5000000000000000]
		])
		if self.modulation:
			pos=np.array([
				[  0.0028042995364173 , 0.9974037352530641 , 0.0077352314347739],
				[  0.9976620421445513 , 0.4989145351814259 , 0.5063342423839297],
				[  0.4945930839448565 , 0.9984163204751321 , 0.5059643310587607],
				[  0.5058732577361121 , 0.4969055205467704 , 0.0045633420079164],
				[  0.5023379578554488 , 0.5015836795248679 , 0.4954366579920836],
				[  0.4971957004635827 , 0.0030944794532296 , 0.9940356689412393],
				[  0.9941267422638879 , 0.5025962647469360 , 0.9936657576160703],
				[  0.0054069160551436 , 0.0010854648185742 , 0.4922647685652261]
			])
		cell=4.8874939611242816*np.eye(3)
		atoms = Atoms('P8',scaled_positions=pos, cell=cell)
		atoms=atoms.repeat([self.latx,self.laty,self.latz])
		atoms.set_pbc([self.xp,self.yp,self.zp])
		return atoms
			
		