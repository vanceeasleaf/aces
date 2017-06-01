from aces.materials  import Material
from ase import Atoms,Atom
from math import pi,sqrt,cos,sin
from aces.materials.graphene import structure as graphene
import numpy as np
class structure(Material):
	def set_parameters(self):
		self.gnrtype='zigzag'
		self.scale=2
		self.line=False
		self.dis=2*self.bond
	def setup(self):
		self.enforceThick=False

	def lmp_structure(self):
		ribbon=graphene(dict(latx=self.latx,laty=self.laty,latz=1,gnrtype=self.gnrtype)).lmp_structure()

		for atom in ribbon:
			atom.position=self.trans(atom.position)
		atoms=Atoms()
		atoms.extend(ribbon)
		atoms.center(vacuum=10)
		
		return atoms

	def trans(self,pos):
		x,y,z=pos
		x1=x
		ds=0.05
		t=0
		if self.line:
			f=self.g
		else:
			f=self.f
		for i in range(int(y/ds)):
			if y<=i*ds:
				break
			r=f(t)
			t+=ds/r
			
		y1=r*cos(t)
		z1=r*sin(t)
		return np.array([x1,y1,z1])

	def f(self,t):
		return self.scale**(t/(2*pi))*self.bond

	def g(self,t):
		return (t/(2*pi)+1)*self.dis

			
		