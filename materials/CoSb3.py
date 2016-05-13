from aces.material import material
from aces.modify import get_unique_atoms
from ase import Atoms,Atom
from math import pi,sqrt
from ase.dft.kpoints import ibz_points
from ase.lattice import bulk
from ase import io
from aces import config
from aces.tools import *
class structure(material):
	def set_parameters(self):
		self.enforceThick=False
		self.latx=1
		self.laty=1
		self.latz=1
		self.elements=['Co','Sb']
		self.poscar="""CoSb3
1.0
        9.0384998322         0.0000000000         0.0000000000
        0.0000000000         9.0384998322         0.0000000000
        0.0000000000         0.0000000000         9.0384998322
   Co   Sb
    8   24
Direct
     0.250000000         0.250000000         0.250000000
     0.750000000         0.750000000         0.750000000
     0.750000000         0.750000000         0.250000000
     0.250000000         0.250000000         0.750000000
     0.250000000         0.750000000         0.750000000
     0.750000000         0.250000000         0.250000000
     0.750000000         0.250000000         0.750000000
     0.250000000         0.750000000         0.250000000
     0.000000000         0.335370004         0.157879993
     0.500000000         0.835370004         0.657880008
     0.000000000         0.664629996         0.842119992
     0.500000000         0.164629996         0.342119992
     0.000000000         0.664629996         0.157879993
     0.500000000         0.164629996         0.657880008
     0.000000000         0.335370004         0.842119992
     0.500000000         0.835370004         0.342119992
     0.157879993         0.000000000         0.335370004
     0.657880008         0.500000000         0.835370004
     0.842119992         0.000000000         0.664629996
     0.342119992         0.500000000         0.164629996
     0.842119992         0.000000000         0.335370004
     0.342119992         0.500000000         0.835370004
     0.157879993         0.000000000         0.664629996
     0.657880008         0.500000000         0.164629996
     0.335370004         0.157879993         0.000000000
     0.835370004         0.657880008         0.500000000
     0.664629996         0.842119992         0.000000000
     0.164629996         0.342119992         0.500000000
     0.335370004         0.842119992         0.000000000
     0.835370004         0.342119992         0.500000000
     0.664629996         0.157879993         0.000000000
     0.164629996         0.657880008         0.500000000
"""
	
	def lmp_structure(self):
		prototype=self.prototype
		atoms=prototype(self.latx,self.laty,self.latz)
		atoms.set_pbc([self.xp,self.yp,self.zp])
		return atoms
		
	
	def prototype(self,latx,laty,latz):
		unit=self.unit()
		col=unit.repeat((latx,laty,latz))
		return col
		
	def unit(self):
		write(self.poscar,"POSCAR_ORI")
		atoms=io.read("POSCAR_ORI")
		return atoms
	

		