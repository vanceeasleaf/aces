from aces.material import material
from ase import Atoms,Atom
from math import pi,sqrt
from aces.materials.graphene import structure as graphene
import numpy as np
class structure(material):
	def set_parameters(self):
		pass
	def setup(self):
		self.xp=self.yp=self.zp=0
		self.enforceThick=False

	def lmp_structure(self):
		atoms=graphene(dict(latx=self.latx,laty=self.laty,latz=self.latz)).lmp_structure()
		center=np.diag(atoms.cell)/2
		for atom in atoms:
			atom.position=self.trans(atom.position,center=center)
		pair=atoms.copy()
		pair.rotate('z',pi/2,center=center)
		pair.rotate('y',pi,center=center)
		pair.translate([0,0,5])
		atoms.extend(pair)
		"""
		atoms.rotate([1,0,1],pi,rotate_cell=True)
		cell=atoms.cell[[2,1,0]]
		cell[1]*=-1
		atoms.set_cell(cell)	
		"""
		self.swap(atoms,1)
		atoms.center(vacuum=10)
		return atoms

	def trans(self,pos,center=[0,0,0]):
		x,y,z=np.array(pos)-center
		dx=0.1
		ds=dx
		y1=y
		x1=0
		sign=x/np.abs(x)
		x=sign*x
		f=self.f
		for i in range(int(x/dx)):
			if x<=i*ds:
				break
			x1+=ds*dx/sqrt((f(x1+dx)-f(x1))**2+dx*dx)
			
		z1=f(x1)
		return np.array([sign*x1,y1,z1])+center

	def f(self,x):
		return .2*x*x
	

			
		