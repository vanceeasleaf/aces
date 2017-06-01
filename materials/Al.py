from aces.materials import Material
from aces.modify import get_unique_atoms
from ase import Atoms,Atom
from math import pi,sqrt
from ase.dft.kpoints import ibz_points
from ase.lattice import bulk
from aces import config
class structure(Material):
	def set_parameters(self):
		self.enforceThick=False
		self.latx=1
		self.laty=1
		self.latz=1
		self.bond=7.653
		self.elements=['Al']
		self.cubic=True

	def setup(self):
		self.bandpoints=ibz_points['fcc']
		self.bandpath=['Gamma','K','X','Gamma','L']
		#self.potential='pair_style	lj/cut %f\n pair_coeff   1 1 %f %f'%(self.cutoff,self.epsilon,self.sigma)
	def lmp_structure(self):
		atoms = bulk('Al', 'fcc', a=self.bond, cubic=self.cubic).repeat([self.latx,self.laty,self.latz])
		atoms.set_pbc([self.xp,self.yp,self.zp])
		atoms.center()
		return atoms
		
	

		