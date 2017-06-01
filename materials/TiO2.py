from aces.materials  import Material
from aces.modify import get_unique_atoms
from ase import Atoms,Atom
from math import pi,sqrt
from ase.dft.kpoints import ibz_points
from ase.lattice.spacegroup import crystal
from aces import config
class structure(Material):
	def set_parameters(self):
		self.enforceThick=False
		self.latx=1
		self.laty=1
		self.latz=1
		self.timestep=self.units.metal.t(.55e-3)
		self.bond=self.units.metal.L(5.430)
		self.elements=['Ti','O']
		self.type='rutile'
	def setup(self):
		#self.bandpoints=ibz_points['fcc']
		#self.bandpath=['Gamma','K','X','Gamma','L']
		#self.potential='pair_style	sw\npair_coeff	* * %s/Si.sw  %s'%(config.lammpspot,' '.join(self.elements))
		pass
	def lmp_structure(self):	
		
		if self.type=='rutile':
			a = 4.584
			c = 2.953
			atoms =crystal(self.elements, basis=[(0, 0, 0), (0.3, 0.3, 0.0)],spacegroup='P4_2/mnm', cellpar=[a, a, c, 90, 90, 90])
		elif self.type=='anatase':
			a=3.7842
			c=9.5146
			atoms=crystal(self.elements, basis=[(0, 0, 0), (0 ,0 ,.2081)],spacegroup=141, cellpar=[a, a, c, 90, 90, 90])
		else:
			raise Exception('unkown TiO2 type!')
		return atoms
		
	

		