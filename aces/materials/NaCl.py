from aces.materials  import Material
from aces.modify import get_unique_atoms
from ase import Atoms,Atom
from math import pi,sqrt
from ase.dft.kpoints import ibz_points
from aces import config
from ase.lattice import bulk
import numpy as np
from aces.tools import *
class structure(Material):
	def set_parameters(self):
		self.cu=False
		pass#['Gamma','Y','T','X','Gamma']

	def setup(self):
		self.forceThick=False
		self.elements=['Na','Cl']
		self.bandpoints=ibz_points['fcc']
		self.bandpath=['Gamma','K','X','Gamma','L']
		self.premitive/=np.array([self.latx,self.laty,self.latz])
		if self.cu:
			self.premitive=np.array([[0,.5,.5],[.5,0,.5],[.5,.5,0]])

	def lmp_structure(self):
		pos=np.array([[0,0,0],[.5,.5,.5]])
		cell=2.8243625205414746*2*np.array([[0,.5,.5],[.5,0,.5],[.5,.5,0]])
		atoms = Atoms('NaCl',scaled_positions=pos, cell=cell)
		if self.cu:
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
			cell=2.8243625205414746*2*np.eye(3)
			atoms = Atoms('Na4Cl4',scaled_positions=pos, cell=cell)
		atoms=atoms.repeat([self.latx,self.laty,self.latz])
		atoms.set_pbc([self.xp,self.yp,self.zp])
		return atoms
			
		