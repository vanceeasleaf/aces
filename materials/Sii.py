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
		self.bond=self.units.metal.L(5.4671121)
		self.elements=['Si']
		self.cubic=True
	def setup(self):
		self.bandpoints=ibz_points['fcc']
		self.bandpath=['Gamma','K','X','Gamma','L']
		self.potential='pair_style	sw\npair_coeff	* * %s/Si.sw  %s'%(config.lammpspot,' '.join(self.elements))
	def lmp_structure(self):
		prototype=self.prototype
		atoms=prototype(self.latx,self.laty,self.latz)
		bond=self.bond
		atoms.set_pbc([self.xp,self.yp,self.zp])
		cell=atoms.cell*bond
		atoms.set_cell(cell,scale_atoms=True)
		return atoms
		
	
	def prototype(self,latx,laty,latz):
		#armchair
		unit=self.unit()
		col=unit.repeat((latx,laty,latz))
		return col
		
	def unit(self):
		cell= [
			[0.0 , 0.5  ,0.5],
	        [0.5  ,0.0  ,0.5],
	        [0.5 , 0.5  ,0.0]
        ]
		atoms=Atoms("Si2",scaled_positions=[[0,0,0],
			[.25,.25,.25]],cell=cell)
		return atoms
	

		