from aces.materials import Material
from aces import default
from ase import Atoms,Atom
class Lead(Material):
	def __init__(self,m):
		self.__dict__=dict(self.__dict__,**default.default)# all the values needed
		self.m=m
		
		self.elements=self.m.elements	
		self.super_setup()
		
	def setup(self):
		self.__dict__=dict(self.__dict__,**self.m.__dict__)
		self.supercell=[2,2,2]
		self.phofc=True
		self.useMini=False

	def lmp_structure(self):
		atoms=Atoms()
		lead1=self.m.lmp_structure()
		self.hatom=len(lead1)
		lead2=lead1.copy()
		lead3=lead1.copy()
		atoms.extend(lead1)
		lead2.translate(lead1.cell[0])
		atoms.extend(lead2)
		lead3.translate(lead1.cell[0]*2)
		atoms.extend(lead3)
		atoms.cell=lead1.cell.copy()
		atoms.cell[0]=lead1.cell[0]*3
		atoms.center()
		return atoms

