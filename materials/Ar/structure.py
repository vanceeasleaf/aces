# encoding : utf8
# C2N hollow 2D structure
from ase import Atoms,Atom
from math import sqrt,pi
from aces import default
from ase.data import atomic_masses,atomic_numbers
from ase.io.vasp import write_vasp
from aces.UnitCell.unitcell import UnitCell
from aces import tools

class structure:
	def __init__(self,home,opt):
		self.home=home
		self.__dict__=dict(self.__dict__,**default.default)# all the values needed
		self.elements=['C','N']
		self.sideLen=3
		self.__dict__=dict(self.__dict__,**opt)
		self.potential='pair_style        tersoff\npair_coeff      * * %s/potentials/BNC.tersoff  %s'%(home,' '.join(self.elements))
		self.masses=""
		self.phontsmasses=""
		i=1
		for a in self.elements:
			num=atomic_numbers[a]
			mass=atomic_masses[num]
			self.masses+="mass %d %f\n"%(i,mass)
			self.phontsmasses+="%s %f 0.0\n"%(a,mass)
			i+=1
		self.dump="dump_modify dump1 element %s"%(' '.join(self.elements))
	def structure(self):
		latx,laty,bond,sideLen=[int(self.latx),int(self.laty),float(self.bond),int(self.sideLen)]
		self.build(latx,laty,bond,sideLen)
		self.write()
		print 'read_data structure'
	def write(self):
		self.atoms.write("structure.xyz")
		write_vasp("POSCAR",self.atoms,sort="True",direct=True,vasp5=True)
		poscar = open("POSCAR")
		unit_cell = UnitCell(poscar)
		unit_cell.num_atom_types=len(self.elements)
		tools.write(unit_cell.output_lammps(),"structure")
import sys		
if __name__=='__main__':
	#structure().build(int(sys.argv[1]),int(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4])).write()
	pass