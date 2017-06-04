from aces.materials  import Material
from ase import Atoms,Atom
from math import pi,sqrt
from ase.dft.kpoints import ibz_points
from ase.lattice import bulk
from aces import config
import numpy as np
class structure(Material):
	def set_parameters(self):
		self.enforceThick=True
		self.thick=1.0
		self.latx=1
		self.laty=1
		self.latz=1
		self.timestep=self.units.metal.t(.55e-3)
		self.elements=['H','He']
		self.fpubeta=1.0
		self.fpug=0.0
		self.creatbonds=1.2
		self.fpua=True
		self.fpub=False
		self.usepre=True
		self.nd=3
		self.ihe=False
		self.fpuk=1.0
	def setup(self):
		if self.usepre:
			self.premitive=np.diag([1.0/self.latx,1.0/self.laty,1.0])
		a='fpu'
		if self.fpua==True:a='fpua'
		if self.fpub:a='fpub'
		#k r0 beta g
		self.potential='neighbor 2.0 nsq\nbond_style	%s\nbond_coeff	* %s 1.0 %s %s\n'%(a,self.fpuk,self.fpubeta,self.fpug)
		#self.potential='bond_style	harmonic\nbond_coeff	* 1.0 1.0\n'
		#self.potential+='\npair_style  none'
	def lmp_structure(self):
		prototype=self.prototype
		atoms=prototype(self.latx,self.laty,self.latz)
		atoms.set_pbc([self.xp,self.yp,self.zp])
		atoms.center()
		if self.ihe:
			symbols=atoms.get_chemical_symbols()
			for i in self.ihe:
				symbols[i]='He'
			atoms.set_chemical_symbols(symbols)
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
	

		