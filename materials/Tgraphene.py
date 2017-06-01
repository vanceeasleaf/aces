from aces.materials  import Material
from aces.modify import get_unique_atoms
from ase import Atoms,Atom
from math import pi,sqrt
from ase.dft.kpoints import ibz_points
class structure(Material):
	def set_parameters(self):
		pass
	def setup(self):
		pass
	def lmp_structure(self):		
		col=self.unitcell(self.laty,self.latx)		
		col.set_pbc([self.xp,self.yp,self.zp])
		atoms=get_unique_atoms(col)
		cell=atoms.cell*self.bond
		atoms.set_cell(cell,scale_atoms=True)
		atoms.center()
		return atoms
		
	
	def unitcell(self,latx,laty):
		unit=Atoms('C2',[(0.5,0.5,0.0),(0.5+sqrt(2)/2,0.5+sqrt(2)/2,0)])
		atoms=Atoms()
		for i in range(4):
			a=unit.copy()
			a.rotate('z',pi/2*i)
			atoms.extend(a)
		atoms.set_cell([2+sqrt(2),2+sqrt(2),10.0])
		col=unit.repeat((latx,laty,1))
		return col
		
			
		