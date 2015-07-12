from aces.material import material
from ase import Atoms,Atom
from math import pi,sqrt,sin,cos
from aces.materials.graphene import structure as graphene
import os
import numpy as np
from aces.tools import *
import aces.config as config
from aces.runners import Runner
class structure(material):
	def set_parameters(self):
		self.scale=5
	def setup(self):
		self.xp=self.yp=self.zp=0
		self.enforceThick=False
		self.useMini=False
	def ori_structure(self):
		atoms=graphene(dict(latx=self.latx,laty=self.laty,latz=1)).lmp_structure()
		lx=self.extent(atoms)[0]
		self.getScale(lx/2)
		self.center_box(atoms)

		low=atoms.positions[:,0].min()+2
		hi=atoms.positions[:,0].max()-2
		self.leftgroup=np.arange(len(atoms),dtype='int')[atoms.positions[:,0]<low]+1
		self.rightgroup=np.arange(len(atoms),dtype='int')[atoms.positions[:,0]>hi]+1
		self.write_graphene(atoms)
		self.posleft=atoms.positions[self.leftgroup[0]-1].copy()+(50,50,50)
		self.posright=atoms.positions[self.rightgroup[0]-1].copy()+(50,50,50)
		self.posleft[0]*=10
		self.posright[0]*=10
		for atom in atoms:
			atom.position=self.trans(atom.position)
		atoms.center(vacuum=50)
		return atoms
	def write_graphene(self,atoms):

		self.atoms=atoms
		mkcd('graphene')
		self.write()
		cd('..')

	def lmp_structure(self):
		if os.path.exists('drag/range'):
			return self.atoms_from_dump('drag/range')
		self.atoms=self.ori_structure()
		drag(self).run()
		return self.atoms_from_dump('drag/range')

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

	def phase(self,x,scale):
		ds=0.1
		dt0=0.01
		t=0
		f=self.px
		for i in range(int(x/ds)):
			if( x<i*ds):break
			vec=f(t+dt0,scale)-f(t,scale)
			t+=ds/np.linalg.norm(vec)*dt0
		return t

	def getScale(self,x):
		# get scale which satisfies phase(x,scale)==2/3*pi,or optimize it
		f=self.phase
		err=1
		ds=0.001
		s=self.scale
		while(err>0.001):
			fd=(f(x,s+ds)-f(x,s))/ds
			dy=pi*7.0/8-f(x,s)
			s+=dy/fd
			err=abs(dy)

		self.scale=s
		debug("scale=%f"%(s))

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
		return self.px(t,self.scale)

	def px(self,t,scale):
		x=sin(t)+2*sin(2*t)
		y=cos(t)-2*cos(2*t)
		z=-sin(3*t)
		return np.array((x,y,z))*scale

class drag(Runner):

	def generate(self):
		m=self.m
		mkdir('drag')
		cd('drag')
		m.write()
		self.getinput()
		passthru(config.mpirun+"  %s "%self.m.cores+config.lammps+" <input  >out.dat")
		cd('..')

	def getinput(self):
		f=open('input', 'w')
		m=self.m
		units,structure,potential,timestep,masses,dumpRate,write_structure,metropolis,useMini,dump=m.units,m.structure,m.potential,m.timestep,m.masses,m.dumpRate,m.write_structure,m.metropolis,m.useMini,m.dump
		print >>f,"units %s"%units
		print >>f,"atom_style atomic"
		print >>f,"boundary s s s"
		print >>f,"dimension 3"
		print >>f,'read_data structure'
		print >>f,potential
		print >>f,"timestep %f"%timestep
		print >>f,masses
		print >>f,"thermo_style custom step pe etotal"
		print >>f,"thermo %d"%dumpRate
		print >>f,"group leftgroup id "+m.toString(m.leftgroup)
		print >>f,"group rightgroup id "+m.toString(m.rightgroup)
		print >>f,"fix dragleft leftgroup drag "+m.toString(m.posleft)+" 1 2"
		print >>f,"fix dragright rightgroup drag "+m.toString(m.posright)+" 1 2"
		print >>f,"reset_timestep 0"
		T=50
		print >>f,"velocity all create %f %d mom yes rot yes dist gaussian"%(T,m.seed)
		print >>f,"fix getEqu  all  nvt temp %f %f %f"%(T,T,m.dtime)

		print >>f,"dump kaka all atom 100 dump.lammpstrj"
		print >>f,"dump_modify  kaka sort id"
		print >>f,"run 200000"
		print >>f,"dump lala all atom 1 range"
		print >>f,"dump_modify  lala sort id"
		print >>f,"run 0"
		f.close()

			
		