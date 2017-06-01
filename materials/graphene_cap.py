from aces.materials  import Material
from ase import Atoms,Atom
from math import pi,sqrt,cos,sin
from aces.materials s.GNT import structure as GNT
from ase.data import extra_molecules
from ase.structure import molecule
import numpy as np
class structure(Material):
	def set_parameters(self):
		pass
	def setup(self):
		self.enforceThick=False
		self.gnrtype='zigzag'
	def lmp_structure(self):
		atoms=GNT(dict(latx=self.latx,laty=5,latz=1,gnrtype=self.gnrtype)).lmp_structure()
		self.center_box(atoms)
		x=atoms.cell[0,0]/2
		c60=molecule('C60',data=extra_molecules.data)
		c60.rotate('z',2*pi/5)
		c60.rotate('y',-pi/4)
		right=c60[c60.positions[:,0]>c60.cell[0,0]/2].copy()
		del c60[c60.positions[:,0]>c60.cell[0,0]/2]
		c60.translate([-x,0,0])
		right.translate([x,0,0])
		atoms.extend(c60)
		atoms.extend(right)
		atoms.center(vacuum=10)	
		return atoms


			
		