from aces.material import material
from aces import default
from ase import Atoms,Atom
class Lead(material):
	def __init__(self,m):
		self.__dict__=dict(self.__dict__,**default.default)# all the values needed
		self.m=m
		self.elements=self.m.elements	
		self.super_setup()
		
	def setup(self):
		self.supercell=[2,2,2]
		self.phofc=True
		self.useMini=False
	def lmp_structure(self):
		atoms=Atoms()
		lead1=self.m.lmp_structure()
		lead2=lead1.copy()
		atoms.extend(lead1)
		lead2.translate(lead1.cell[0])
		atoms.extend(lead2)
		atoms.cell=lead1.cell.copy()
		atoms.cell[0]=lead1.cell[0]+lead2.cell[0]
		atoms.center()
		return atoms

