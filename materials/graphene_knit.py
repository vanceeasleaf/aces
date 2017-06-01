from aces.materials  import Material
from ase import Atoms,Atom
from math import pi,sqrt,cos,sin
from aces.materials s.graphene import structure as graphene
import numpy as np
from aces.tools import *
class structure(Material):
	def set_parameters(self):
		self.amp=1.7
	def setup(self):
		self.enforceThick=False

	def lmp_structure(self):
		ribbon=graphene(dict(latx=20,laty=3,latz=1,gnrtype='zigzag')).lmp_structure()
		
		self.length=ribbon.cell[0,0]
		self.center_box(ribbon)
		self.w=.2
		self.getW(self.length/2)
		newlx=2*self.phase(self.length/2,self.w)
		for atom in ribbon:
			atom.position=self.trans(atom.position)
		atoms=ribbon
		#atoms.center(vacuum=10,axis=[1,2])
		#atoms.cell[0,0]=newlx
		#atoms.center(axis=[0])
		newatoms=Atoms()
		for i in range(4):
			atoms.rotate('z',i*pi/2,center=[0,-newlx/4,0])
			newatoms.extend(atoms.copy())
		newatoms.cell[0,0]=newatoms.cell[1,1]=newlx
		newatoms.center()
		atoms=newatoms#.repeat([self.latx,self.laty,1])
		return atoms

	def getW(self,x):
		f=lambda s:self.phase(x,s)*s
		err=1
		ds=0.001
		s=self.w
		n=0
		while(err>0.001 and n<1000):
			fd=(f(s+ds)-f(s))/ds
			dy=pi-f(s)
			s+=dy/fd
			err=abs(dy)
			n+=1
		self.w=s
		debug("w=%f,f(s)=%f"%(s,f(s)))

	def trans(self,pos):
		x,y,z=pos
		y1=y
		t=self.phase(x,self.w)
		f=self.f
		xc,yc,zc=f(t)
		vec=self.tanvec(f,t)
		x1,y1,z1=np.array([xc,yc,zc])+y*self.yvec(f,t)+z*self.normvec(f,t)
		return np.array([x1,y1,z1])

	def f(self,t):
		return np.array([t,0,self.amp*sin(self.w*t)])

	def curve(self,t,w):
		return np.array([t,0,self.amp*sin(w*t)])

	def phase(self,x,w):
		ds=0.1
		dt0=0.01
		t=0
		f=self.curve
		sign=x/abs(x)
		x=abs(x)
		dt0*=sign
		for i in range(int(x/ds)):
			if( x<i*ds):break
			vec=f(t+dt0,w)-f(t,w)
			t+=ds/np.linalg.norm(vec)*dt0
		return t		
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