from aces.material import material
from aces import default
from ase import Atoms,Atom
class Device(material):
	def __init__(self,m,m1,m2):
		self.__dict__=dict(self.__dict__,**default.default)# all the values needed
		self.m=m
		self.m1=m1
		self.m2=m2
		self.elements=self.m.elements	
		self.super_setup()
		
	def setup(self):
		self.supercell=[1,1,1]
		self.phofc=True
		self.useMini=False
	def lmp_structure(self):
		atoms=Atoms()
		lead1=self.m1.lmp_structure()
		center=self.m.lmp_structure()
		lead2=self.m2.lmp_structure()
		atoms.extend(lead1)
		center.translate(lead1.cell[0])
		atoms.extend(center)
		lead2.translate(lead1.cell[0]+center.cell[0])
		atoms.extend(lead2)
		atoms.set_pbc([0,0,0])
		atoms.cell[0]=lead1.cell[0]+center.cell[0]+lead2.cell[0]
		atoms.center(10.0,axis=[1,2])
		x=atoms.positions[:,0]

		return atoms

