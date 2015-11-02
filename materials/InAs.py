from aces.material import material
from aces.modify import get_unique_atoms
from ase import Atoms,Atom
from math import pi,sqrt
from ase.dft.kpoints import ibz_points
class structure(material):
	def set_parameters(self):
		self.bond=6.059141
	def setup(self):
		self.elements=["In","As"]
		self.bandpoints=ibz_points['fcc']
		self.bandpath=['Gamma','K','X','Gamma','L']
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
		atoms=Atoms("InAs",scaled_positions=[[0,0,0],
			[.25,.25,.25]],cell=cell)
		return atoms
			
		