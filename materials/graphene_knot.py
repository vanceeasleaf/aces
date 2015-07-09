from aces.material import material
from ase import Atoms,Atom
from math import pi,sqrt,sin,cos
from aces.materials.graphene import structure as graphene
import numpy as np
class structure(material):
	def set_parameters(self):
		self.scale=5
	def setup(self):
		self.xp=self.yp=self.zp=0
		self.enforceThick=False

	def lmp_structure(self):
		atoms=graphene(dict(latx=self.latx,laty=self.laty,latz=self.latz)).lmp_structure()
		self.center_box(atoms)
		for atom in atoms:
			atom.position=self.trans(atom.position)
		atoms.center(vacuum=10)
		return atoms

	def trans(self,pos):
		x,y,z=pos
		ds=0.1
		dt0=0.01
		t=0
		f=self.knot
		sign=x/abs(x)
		x=abs(x)
		dt0*=sign
		for i in range(int(x/ds)):
			if( x<i*ds):break
			vec=f(t+dt0)-f(t)
			t+=ds/np.linalg.norm(vec)*dt0
		xc,yc,zc=f(t)
		vec=self.tanvec(f,t)
		x1,y1,z1=np.array([xc,yc,zc])+y*self.yvec(f,t)+z*self.normvec(f,t)
		return np.array([x1,y1,z1])

	def tanvec(self,f,t):
		dt0=0.01
		vec=(f(t+dt0)-f(t))/dt0
		vec/=np.linalg.norm(vec)
		return vec

	def normvec(self,f,t):
		dt0=0.01
		v=self.tanvec
		n=(v(f,t+dt0)-v(f,t))/dt0
		n/=np.linalg.norm(n)
		return n

	def yvec(self,f,t):
		v=self.tanvec(f,t)
		n=self.normvec(f,t)
		return np.cross(n,v)

	def knot(self,t):
		x=sin(t)+2*sin(2*t)
		y=cos(t)-2*cos(2*t)
		z=-sin(3*t)
		return np.array((x,y,z))*self.scale

	

			
		