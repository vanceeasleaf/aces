from aces.material import material
from ase import Atoms,Atom
from math import pi,sqrt,cos,sin
from aces.materials.graphene import structure as graphene
import numpy as np
from aces.runners import Runner
from aces.runners.phonopy import RotateVector
from aces.tools import *
import aces.config as config
import os
class structure(material):
	def set_parameters(self):
		self.angle=pi/4
		self.twist=0.0
	def setup(self):
		self.xp=self.yp=self.zp=0
		self.enforceThick=False
		self.useMini=False
	def ori_structure(self):
		atoms=graphene(dict(latx=self.latx,laty=self.laty,latz=self.latz)).lmp_structure()
		self.center_box(atoms)
		for atom in atoms:
			atom.position=self.trans(atom.position)
		pair=atoms.copy()
		pair.rotate('z',pi/2)
		pair.rotate('y',pi)
		pair.translate([0,0,10])
		atoms.extend(pair)
		self.swap(atoms,1)
		dx=7
		xlo=atoms.positions[:,0].min()+dx
		xhi=atoms.positions[:,0].max()-dx
		ylo=atoms.positions[:,1].min()+dx
		yhi=atoms.positions[:,1].max()-dx
		zlo=atoms.positions[:,2].min()+dx
		zhi=atoms.positions[:,2].max()-dx
		self.x1z1=np.arange(len(atoms),dtype='int')[np.logical_and(atoms.positions[:,0]>xhi , atoms.positions[:,2]>zhi)]+1
		self.x1z0=np.arange(len(atoms),dtype='int')[np.logical_and(atoms.positions[:,0]>xhi , atoms.positions[:,2]<zlo)]+1
		self.x0y1=np.arange(len(atoms),dtype='int')[np.logical_and(atoms.positions[:,0]<xlo , atoms.positions[:,1]>yhi)]+1
		self.x0y0=np.arange(len(atoms),dtype='int')[np.logical_and(atoms.positions[:,0]<xlo , atoms.positions[:,1]<ylo)]+1
		g=self.angle/2.0
		self.posx1z1=np.array((50,50,50))+500*np.array([cos(g),0.0,sin(g)])
		self.posx1z0=np.array((50,50,50))+500*np.array([cos(g),0.0,-sin(g)])
		self.posx0y1=np.array((50,50,50))+500*np.array([-cos(g),sin(g),0.0])
		self.posx0y0=np.array((50,50,50))+500*np.array([-cos(g),-sin(g),0.0])
		atoms.center(vacuum=50)
		return atoms
	def lmp_structure(self):
		if os.path.exists('drag/range'):
			return self.atoms_from_dump('drag/range')
		self.atoms=self.ori_structure()
		drag(self).run()
		return self.atoms_from_dump('drag/range')
	def trans(self,pos):
		x,y,z=np.array(pos)
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
		return np.array([sign*x1,y1,z1])

	def f(self,x):
		return .1*x*x
	
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
		print >>f,"group x1z0 id "+m.toString(m.x1z0)
		print >>f,"group x1z1 id "+m.toString(m.x1z1)
		print >>f,"group x0y0 id "+m.toString(m.x0y0)
		print >>f,"group x0y1 id "+m.toString(m.x0y1)

		print >>f,"reset_timestep 0"
		T=50
		print >>f,"velocity all create %f %d mom yes rot yes dist gaussian"%(T,m.seed)
		print >>f,"fix getEqu  all  nvt temp %f %f %f"%(T,T,m.dtime)
		print >>f,"fix 1 all recenter 50 50 50 units box"
		print >>f,"dump kaka all atom 100 dump.lammpstrj"
		print >>f,"dump_modify  kaka sort id"
		print >>f,"fix dragx0y0 x0y0 drag "+m.toString(m.posx0y0)+" .5 2"
		print >>f,"fix dragx0y1 x0y1 drag "+m.toString(m.posx0y1)+" .5 2"
		for i in range(20):
			g=m.angle/2.0
			m.posx1z1=np.array((50,50,50))+500*RotateVector(np.array([cos(g),0.0,sin(g)]),[1,0,0],m.twist/20.0*i)
			m.posx1z0=np.array((50,50,50))+500*RotateVector(np.array([cos(g),0.0,-sin(g)]),[1,0,0],m.twist/20.0*i)
			print >>f,"fix dragx1z0 x1z0 drag "+m.toString(m.posx1z0)+" .5 2"
			print >>f,"fix dragx1z1 x1z1 drag "+m.toString(m.posx1z1)+" .5 2"
			print >>f,"run 10000"
			print >>f,"unfix dragx1z0"
			print >>f,"unfix dragx1z1"
		print >>f,"dump lala all atom 1 range"
		print >>f,"dump_modify  lala sort id"
		print >>f,"run 0"
		f.close()

			
		