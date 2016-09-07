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
		pass#['Gamma','Y','T','X','Gamma']

	def setup(self):
		self.forceThick=False
		self.elements=['Si','C']
		self.bandpoints=ibz_points['fcc']
		self.bandpath=['Gamma','K','X','Gamma','L']
		self.premitive/=np.array([self.latx,self.laty,self.latz])

	def lmp_structure(self):
		pos=np.array([[0,0,0],[.25,.25,.25]])
		cell=8.237*0.5291772083*np.array([[0,.5,.5],[.5,0,.5],[.5,.5,0]])
		atoms = Atoms('SiC',scaled_positions=pos, cell=cell)
		atoms=atoms.repeat([self.latx,self.laty,self.latz])
		atoms.set_pbc([self.xp,self.yp,self.zp])
		return atoms
			
		