from aces.material import material
from aces.modify import get_unique_atoms
from ase import Atoms,Atom
from math import pi,sqrt
from ase.dft.kpoints import ibz_points
from aces import config
import numpy as np
from aces.tools import *
class structure(material):
	def set_parameters(self):
		self.gnrtype='armchair'
		self.thick=20.0
		self.choice=0

	def setup(self):
		es=[
			['Sn','S' ],
			['Sn','Se'],
			['Ge','S' ],
			['Ge','Se'],
			['Sn','Te']
		]
		self.elements=es[self.choice]
		self.bandpoints=ibz_points['orthorhombic']
		self.bandpoints['T']=self.bandpoints['S']
		self.bandpath=['Gamma','Y','T','X','Gamma']
		self.premitive/=np.array([self.latx,self.laty,self.latz])

	def lmp_structure(self):
		prototype=self.prototype
		if self.gnrtype=='armchair':
			col=prototype(self.latx,self.laty)		
		elif self.gnrtype=='zigzag':			
			col=prototype(self.laty,self.latx)		
			self.swap(col,2)
		else: raise Exception('Unkown gnr type!')	
		col.center()	
		atoms=get_unique_atoms(col)
		col.set_pbc([self.xp,self.yp,self.zp])
		atoms.center()
		return atoms
		
	
	
	def prototype(self,latx,laty):
		#armchair
		atoms=self.unit()
		col=atoms.repeat((latx,laty,1))
		return col
		
	def unit(self):
		labels=toString(np.repeat(self.elements,2),'')
		pos=np.array([  0.0890026531505678  ,0.25 ,0.4473684554803318,
		  0.5890026531501628 , 0.75  ,0.5526315445196681,
  0.4109973468498373 , 0.75 , 0.4473684554803318,
  0.9109973468498372 , 0.25, 0.5526315445196681]).reshape([4,3])
		cell=np.diag([4.61,3.299,self.thick])
		atoms=Atoms(labels,scaled_positions=pos,cell=cell)
		#atoms.set_cell(cell)
		atoms.center()
		#atoms.set_cell([4.35,4.02,self.thick])
		return atoms
			
		