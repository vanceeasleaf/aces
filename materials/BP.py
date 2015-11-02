from aces.material import material
from aces.modify import get_unique_atoms
from ase import Atoms,Atom
from math import pi,sqrt
from ase.dft.kpoints import ibz_points
from aces import config
class structure(material):
	def set_parameters(self):
		self.gnrtype='armchair'
		self.ax=4.6059999466
		self.ay=3.3002998829
		self.thick=20.0
		self.bond=2.22136
		self.bond1=2.25747
		self.pot=False
	def setup(self):
		self.elements=["P"]
		self.bandpoints=ibz_points['orthorhombic']
		self.bandpath=['Gamma','Y','S','X','Gamma']
		self.t=sqrt(self.bond**2-(self.ay/2.0)**2)
		self.ox=(self.ax-2.0*self.t)/2.0
		self.oz=sqrt(self.bond1**2-(self.ox)**2)
		if self.pot==1:
			print "BP potential chosen:P.sw"
			self.potential='pair_style	sw\npair_coeff	* * %s/P.sw  %s'%(config.lammpspot,' '.join(self.elements))

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
		t=self.t
		atoms=Atoms('P4',
			[(0.0,0.0,0.0),
			(t,self.ay*0.5,0.0),
			(t+self.ox,self.ay*0.5,self.oz),
			(self.ax-self.ox,0.0,self.oz)
			])
		atoms.set_cell([self.ax,self.ay,self.thick])
		return atoms
			
		