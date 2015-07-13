from aces.material import material
from ase import Atoms,Atom
from math import pi,sqrt,cos,sin
from aces.materials.GNT import structure as GNT
import numpy as np
class structure(material):
	def set_parameters(self):
		self.gnrtype='zigzag'
		self.angle=2*pi
		self.phi=2*pi
	def setup(self):
		self.enforceThick=False

	def lmp_structure(self):
		atoms=GNT(dict(latx=self.latx,laty=self.laty,latz=1,gnrtype=self.gnrtype,phi=self.angle)).lmp_structure()
		r=atoms.cell[0,0]/self.phi
		gnt=GNT()
		self.center_box(atoms)
		atoms.rotate('z',pi/2)
		for atom in atoms:
			atom.position=gnt.trans(atom.position,r=r)
		#atoms.rotate('y',pi/2)
		atoms.center(vacuum=10)	
		return atoms


			
		