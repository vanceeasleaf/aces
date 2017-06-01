from aces.materials  import Material
from ase import Atoms,Atom
from math import pi,sqrt,cos,sin
from aces.materials s.graphene_twist import structure as Twist
from aces.materials s.GNT import structure as GNT
import numpy as np
class structure(Material):
	def set_parameters(self):
		self.gnrtype='zigzag'
	def setup(self):
		self.enforceThick=False

	def lmp_structure(self):
		atoms=Twist(dict(latx=self.latx,laty=self.laty,latz=1,gnrtype=self.gnrtype,twist=pi)).lmp_structure()
		self.writeatoms(atoms,'twist')
		r=atoms.cell[0,0]/2/pi
		gnt=GNT()
		self.center_box(atoms)
		atoms.rotate('z',pi/2)
		for atom in atoms:
			atom.position=gnt.trans(atom.position,r=r)
		#atoms.rotate('y',pi/2)
		atoms.center(vacuum=10)	
		return atoms


			
		