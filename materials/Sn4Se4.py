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
		pass#['Gamma','Y','T','X','Gamma']

	def setup(self):
		self.forceThick=False
		self.elements=['Se','Sn']
		self.bandpoints=ibz_points['orthorhombic']
		self.bandpoints['T']=self.bandpoints['S']
		self.bandpath=["Gamma","X","Z","Gamma","Y","S","R","Gamma"]
		self.premitive/=np.array([self.latx,self.laty,self.latz])

	def lmp_structure(self):
		prototype=self.prototype
		col=prototype(self.latx,self.laty,self.latz)		
		atoms=get_unique_atoms(col)
		col.set_pbc([self.xp,self.yp,self.zp])
		return atoms
		
	
	
	def prototype(self,latx,laty,latz):
		atoms=self.unit()
		col=atoms.repeat((latx,laty,latz))
		return col
		
	def unit(self):
		pos=np.array([
				[0.75000000000000,   0.48503746453728,   0.85560968541486],   
				[0.25000000000000,   0.98503746453728,   0.64439031458514],   
				[0.25000000000000,   0.51496253546272,   0.14439031458514],   
				[0.75000000000000,   0.01496253546272,   0.35560968541486],   
				[0.75000000000000,   0.09621535246053,   0.12506317936786],   
				[0.25000000000000,   0.59621535246053,   0.37493682063214],   
				[0.25000000000000,   0.90378464753947,   0.87493682063214],   
				[0.75000000000000,   0.40378464753947,   0.62506317936786]
   ])
		cell=np.diag([4.27621131928033,4.54139940515969,11.98738752160755])
		atoms=Atoms('Se4Sn4',scaled_positions=pos,cell=cell)
		#atoms.set_cell(cell)
		#atoms.center()
		#atoms.set_cell([4.35,4.02,self.thick])
		return atoms
			
		