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
		self.latx=1
		self.laty=1
		self.latz=1
		self.timestep=self.units.metal.t(.55e-3)
		self.bond=self.units.metal.L(5.430)
		self.elements=['C']
		self.cubic=True
		self.al=False
	def setup(self):
		self.bandpoints=ibz_points['fcc']
		self.bandpath=['Gamma','K','X','Gamma','L']
		self.potential='pair_style	sw\npair_coeff	* * %s/Si.sw  %s'%(config.lammpspot,' '.join(self.elements))
	def lmp_structure(self):
		atoms = bulk('C', 'diamond', a=self.bond, cubic=self.cubic)
		atoms=atoms.repeat([self.latx,self.laty,self.latz])
		atoms.set_pbc([self.xp,self.yp,self.zp])
		
		if not self.al:
			atoms.center()
		return atoms
		
	

		