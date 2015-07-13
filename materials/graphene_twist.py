from aces.material import material
from ase import Atoms,Atom
from math import pi,sqrt,cos,sin
from aces.materials.graphene import structure as graphene
import numpy as np
class structure(material):
	def set_parameters(self):
		self.gnrtype='zigzag'
		self.twist=pi
	def setup(self):
		self.enforceThick=False

	def lmp_structure(self):
		ribbon=graphene(dict(latx=self.latx,laty=self.laty,latz=1,gnrtype=self.gnrtype)).lmp_structure()
		self.writeatoms(ribbon,'ribbon')
		self.length=ribbon.cell[0,0]
		self.center_box(ribbon)
		for atom in ribbon:
			atom.position=self.trans(atom.position)
		atoms=ribbon
		atoms.center(vacuum=10,axis=[1,2])
		atoms.center(axis=[0])
		return atoms

	def trans(self,pos):
		x,y,z=pos
		x1=x
		t=self.twist/self.length*x
		y1=y*cos(t)-z*sin(t)
		z1=y*sin(t)+z*cos(t)
		return np.array([x1,y1,z1])


			
		