from aces.materials  import Material
from ase import io
from ase import Atoms,Atom
import numpy as np
from aces.tools import *
class structure(Material):
	def setup(self):
		self.forceThick=False		
		self.POSCAR=self.getPOSCAR()

		self.setElements()
		self.premitive/=np.array([self.latx,self.laty,self.latz])
		self.csetup()
	def csetup(self):
		pass
	def setElements(self):
		self.elements=self.POSCAR.split('\n')[5].strip().split()
	def lmp_structure(self):
		write(self.POSCAR,"poscar_aces_structure")
		atoms=io.read("poscar_aces_structure",format="vasp")
		atoms=atoms.repeat([self.latx,self.laty,self.latz])
		atoms.set_pbc([self.xp,self.yp,self.zp])
		return atoms
			
		