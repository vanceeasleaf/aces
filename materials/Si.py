from aces.material import material
from aces.modify import get_unique_atoms
from ase import Atoms,Atom
from math import pi,sqrt
from ase.dft.kpoints import ibz_points
from ase.lattice import bulk
from aces import config
class structure(material):
	def set_parameters(self):
		self.enforceThick=False
		self.latx=288
		self.laty=4
		self.latz=4
		self.timestep=self.units.metal.t(.55e-3)
		self.bond=self.units.metal.L(5.430)
		self.elements=['Si']
		self.cubic=True
	def setup(self):
		self.bandpoints=ibz_points['fcc']
		self.bandpath=['Gamma','K','X','Gamma','L']
		self.potential='pair_style	tersoff\npair_coeff	* * %s/SiC_1994.tersoff  %s'%(config.lammpspot,' '.join(self.elements))
	def lmp_structure(self):
		atoms = bulk('Si', 'diamond', a=self.bond, cubic=self.cubic)
		atoms.set_pbc([self.xp,self.yp,self.zp])
		atoms.center()
		return atoms
		
	

		