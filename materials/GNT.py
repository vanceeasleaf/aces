from aces.material import material
from ase import Atoms,Atom
from math import pi,sqrt,cos,sin
from aces.materials.graphene import structure as graphene
import numpy as np
class structure(material):
	def set_parameters(self):
		self.type='zigzag'
		self.phi=2*pi
	def setup(self):
		self.enforceThick=False

	def lmp_structure(self):
		atoms=graphene(dict(latx=self.latx,laty=self.laty,latz=self.latz,gnrtype=self.type)).lmp_structure()
		center=np.diag(atoms.cell)/2
		self.radius=atoms.cell[1][1]/self.phi
		for atom in atoms:
			atom.position=self.trans(atom.position,center=center,r=self.radius)

		atoms.center(vacuum=10,axis=[1,2])
		atoms.center(vacuum=0,axis=[0])
		return atoms

	def trans(self,pos,center=[0,0,0],r=1):
		x,y,z=np.array(pos)-center
		x1=x
		t=y/r
		y1=sin(t)*r
		z1=r-cos(t)*r
		return np.array([x1,y1,z1])+center

	

			
		