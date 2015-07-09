from aces.material import material
from aces.modify import get_unique_atoms
from ase import Atoms,Atom
from math import pi,sqrt
from ase.dft.kpoints import ibz_points
class structure(material):
	def set_parameters(self):
		self.centerN=False
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
		unit=Atoms('C2N',[(1/2*sqrt(3),0.5,0.0),
			(1/2*sqrt(3),-0.5,0),
			(1/2*sqrt(3),1.5,0)])
		atoms=Atoms()
		for i in range(3):
			a=unit.copy()
			a.rotate('z',pi*2/3*i)
			atoms.extend(a)
		if self.centerN:
			atoms.append(Atom('N',(0,0,0)))
		atoms.set_cell([6*sqrt(3)/2,6,10.0])
		col=unit.repeat((latx,laty,1))
		return col
		
			
		