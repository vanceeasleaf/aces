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
		self.elements=['Si','S']
		self.cubic=False
		self.CUBIC="""cubic sulfied\silicene
	 1.00000000000000		 
		 6.6442279437829646		0.0000000000000000		0.0000000000000000
		 0.0000000000000000		3.9854738543450634		0.0000000000000000
		 0.0000000000000000		0.0000000000000000		9.4373464922296861
	 Si	 S 
		 4		 4
Direct
	0.0000000000000000	0.0000000000000000	0.5388405482746847
	0.5000000000000000	0.0000000000000000	0.5388405482746847
	0.0000000000000000	0.5000000000000000	0.6778181711064306
	0.5000000000000000	0.5000000000000000	0.6778181711064306
	0.2500000000000000	0.9547368054194507	0.3925572426466033
	0.7500000000000000	0.0452631945805493	0.3925572426466033
	0.2500000000000000	0.5448843497126532	0.8241070379723041
	0.7500000000000000	0.4551156502873468	0.8241070379723041
		"""
		self.NONCUBIC="""sulfied\silicene												
	 1.00000000000000		 
		 4.5575267218017510		0.0016234095507666		0.0000000000000000
		-2.7951936058658609		3.5997145943922164		0.0000000000000000
		 0.0000000000000000		0.0000000000000000		8.4548765236271368
	 Si	 S 
		 2		 2
Direct
	0.1703001548270125	0.8296997861729879	0.5000000000000000
	0.8296998451729876	0.1703002138270122	0.5000000000000000
	0.0000000000000000	0.0000000000000000	0.6968304092241721
	0.0000000000000000	0.0000000000000000	0.3031695907758277
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
		if self.cubic:
			write(self.CUBIC,"POSCAR_ORI")
		else:
			write(self.NONCUBIC,"POSCAR_ORI")
		atoms=io.read("POSCAR_ORI")
		return atoms
	

		