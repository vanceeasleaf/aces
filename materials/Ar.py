from aces.materials  import Material
from aces.modify import get_unique_atoms
from ase import Atoms,Atom
from math import pi,sqrt
from ase.dft.kpoints import ibz_points
from ase.lattice import bulk
from aces import config
class structure(Material):
	def set_parameters(self):
		self.enforceThick=False
		self.latx=4
		self.laty=4
		self.latz=4
		self.timestep=self.units.metal.t(.55e-3)
		self.bond=self.units.metal.L(5.260)
		self.elements=['Ar']
		self.cubic=True
		self.cutoff=self.units.lj.L(2.8)
		self.epsilon=self.units.lj.E(1)
		self.sigma=self.units.lj.L(1)

	def setup(self):
		self.bandpoints=ibz_points['fcc']
		self.bandpath=['Gamma','K','X','Gamma','L']
		self.potential='pair_style	lj/cut %f\n pair_coeff   1 1 %f %f'%(self.cutoff,self.epsilon,self.sigma)
	def lmp_structure(self):
		atoms = bulk('Ar', 'fcc', a=self.bond, cubic=self.cubic).repeat([self.latx,self.laty,self.latz])
		atoms.set_pbc([self.xp,self.yp,self.zp])
		atoms.center()
		return atoms
		
	

		