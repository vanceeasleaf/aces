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
		self.angle=pi/4
	def setup(self):
		self.enforceThick=False
		self.gnrx=int(self.latx/cos(self.angle))+2
	def lmp_structure(self):
		ribbon=graphene(dict(latx=self.gnrx,laty=self.gnry,latz=1,gnrtype=self.gnrtype)).lmp_structure()
		gnt=GNT(dict(latx=self.latx,laty=self.laty,latz=1,type=self.gnttype))
		atoms=gnt.lmp_structure()
		self.center_box(atoms)
		
		self.center_box(ribbon)
		ribbon.rotate('z',self.angle)
		r=gnt.radius+self.bond
		for atom in ribbon:
			atom.position=gnt.trans(atom.position,r=r)
		max=atoms.positions[:,0].max()
		min=atoms.positions[:,0].min()
		del ribbon[ribbon.positions[:,0]>max]
		del ribbon[ribbon.positions[:,0]<min]
		ribbon.translate([0,0,-r])
		atoms.extend(ribbon)
		atoms.center(vacuum=10)
		
		return atoms

	

			
		