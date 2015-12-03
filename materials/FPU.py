from aces.material import material
from aces.modify import get_unique_atoms
from ase import Atoms,Atom
from math import pi,sqrt
from ase.dft.kpoints import ibz_points
from ase.lattice import bulk
from aces import config
import numpy as np
class structure(material):
	def set_parameters(self):
		self.enforceThick=True
		self.thick=1.0
		self.latx=1
		self.laty=1
		self.latz=1
		self.timestep=self.units.metal.t(.55e-3)
		self.elements=['H']
		self.fpubeta=1.0
		self.fpug=0.0
		self.creatbonds=1.2
		self.fpua=False
		self.fpub=False
		self.usepre=True
		self.nd=3
	def setup(self):
		if self.usepre:
			self.premitive=np.diag([1.0/self.latx,1.0/self.laty,1.0])
		a='fpu'
		if self.fpua==True:a='fpua'
		if self.fpub:a='fpub'
		#k r0 beta g
		self.potential='bond_style	%s\nbond_coeff	* 1.0 1.0 %s %s\n'%(a,self.fpubeta,self.fpug)
		#self.potential='bond_style	harmonic\nbond_coeff	* 1.0 1.0\n'
		#self.potential+='\npair_style  none'
	def lmp_structure(self):
		prototype=self.prototype
		atoms=prototype(self.latx,self.laty,self.latz)
		atoms.set_pbc([self.xp,self.yp,self.zp])
		atoms.center()
		return atoms
		
	
	def prototype(self,latx,laty,latz):
		#armchair
		unit=self.unit()
		col=unit.repeat((latx,laty,latz))
		return col
		
	def unit(self):
		nz=10
		ny=1
		if not self.enforceThick:
			nz=1
		if self.dimension==1:
			ny=10
			nz=10
		elif self.dimension==2:
			ny=1
			nz=10
		cell=np.diag([1,ny,nz]) 
		atoms=Atoms("H",scaled_positions=[[0,0,0]],cell=cell)
		return atoms
	

		