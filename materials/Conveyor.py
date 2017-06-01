from aces.materials  import Material
from ase import Atoms,Atom
from math import pi,sqrt,cos,sin
from aces.materials s.graphene import structure as graphene
from aces.materials s.GNT import structure as GNT
import numpy as np
class structure(Material):
	def set_parameters(self):
		self.gnrtype='zigzag'
		self.gnttype='zigzag'
		self.gnry=2
		self.gnrx=20
		self.gntx=3
		self.gnty=4

	def setup(self):
		self.enforceThick=False

	def lmp_structure(self):
		ribbon=graphene(dict(latx=self.gnrx,laty=self.gnry,latz=1,gnrtype=self.gnrtype)).lmp_structure()
		length=ribbon.cell[0][0]
		gnt=GNT(dict(latx=self.gntx,laty=self.gnty,latz=1,type=self.gnttype))
		atoms=gnt.lmp_structure()
		self.center_box(atoms)
		self.center_box(ribbon)
		atoms.rotate('z',pi/2)
		left=atoms.copy()

		r=gnt.radius+self.bond
		d=(length-2*pi*r)/2
		left.translate([-d/2,0,0])
		atoms.translate([d/2,0,0])
		for atom in ribbon:
			atom.position=self.trans(atom.position,r=r,d=d)
		atoms.extend(ribbon)
		atoms.extend(left)
		atoms.center(vacuum=10)
		
		return atoms

	def trans(self,pos,r,d):
		x,y,z=pos
		y1=y
		if abs(x)<d/2:
			x1=x			
			z1=z+r
		elif abs(x)<d/2+pi*r:
			t=(abs(x)-d/2)/r
			x1=d/2+sin(t)*r
			x1*=x/abs(x)
			z1=cos(t)*r
		else:
			xs=d/2+pi*r
			x1=d/2-(abs(x)-xs)
			x1*=x/abs(x)
			z1=z-r
		return x1,y1,z1
			
		