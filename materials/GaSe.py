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
		pass

	def setup(self):
		self.elements=['Ga','Se']
		self.bandpoints=ibz_points['hexagonal']
		self.bandpath=['Gamma','M','K','L','Gamma']
		self.premitive/=np.array([self.latx,self.laty,self.latz])

	def lmp_structure(self):
		pos=np.array([[1.9091191292,2.2044608863,	14.5257967802],
					  [1.9091191292,2.2044608863,	12.0542563707],
					  [0.0000000000,1.1022304431,	10.8813264123],
					  [0.0000000000,1.1022304431,	15.6992399940]])


		cell=np.array([
			[3.8182382584	,0.0000000000	,0.0000000000	],
		    [-1.9091191292	,3.3066913294	,0.0000000000	],
		    [0.0000000000	,0.0000000000	,26.5805664063	]])
		atoms = Atoms('Ga2Se2',positions=pos, cell=cell)
		atoms=atoms.repeat([self.latx,self.laty,self.latz])
		atoms.set_pbc([self.xp,self.yp,self.zp])
		return atoms
			
		