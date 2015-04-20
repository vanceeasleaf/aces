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
	def build(self,latx,laty,bond,sideLen):
		self.sideLen=sideLen
		self.ax=(sideLen+1)*3*sqrt(3)/2
		self.ay=(sideLen+1)*1.5*2
		unitCell=self.quasiCell(latx,laty)
		unitCell.set_cell([latx*self.ax,laty*self.ay,100])
		unitCell.set_cell([latx*self.ax*bond,laty*self.ay*bond,100],scale_atoms=True)
		unitCell.set_pbc([1,1,1])
		unitCell.center()
		sorted(unitCell,key=lambda atom : atom.position[0])
		self.atoms=unitCell
	def phontsAtoms(self):
		cell=self.atoms.get_cell()
		content="cell %f %f %f\n"%(cell[0][0],cell[1][1],cell[2][2])
		content+="natoms %d\n"%(len(self.atoms))
		content+="fractional\n"
		pos=self.atoms.get_scaled_positions()
		for i,atom in enumerate(self.atoms):
			content+="%s %s\n"%(atom.symbol,' '.join(["%s"%x for x in pos[i]]))
		return content
	def write(self):
		self.atoms.write("structure.xyz")
		write_vasp("POSCAR",self.atoms,sort="True",direct=True,vasp5=True)
		poscar = open("POSCAR")
		unit_cell = UnitCell(poscar)
		unit_cell.num_atom_types=len(self.elements)
		tools.write(unit_cell.output_lammps(),"structure")
		#tools.write(self.phontsAtoms(),"phonts.data")

		#data("structure").write("new")
		
	def quasiCell(self,col,row,type=0):
		atoms=self.getUnitCell()
		unitCell=Atoms()
		for i in range(col):
			for j in range(row):
				x=i*self.ax;
				y=(j+(type-0.5)*(i%2))*self.ay;
				newAtoms=atoms.copy()
				newAtoms.translate([x,y,0])
				unitCell.extend(newAtoms)
		return unitCell
		
	def getUnitCell(self):
		sideLen=self.sideLen	
		edge=Atoms()
		nAtom=2*sideLen+1
		for i in range(nAtom):
			if i%2==0:
				label='C'
				y=0.5
			else:
				label='N'
				if len(self.elements)==1:label='C'
				y=1.0
			x=(-sideLen+i)*0.5*sqrt(3)
			y-=(sideLen+1)*1.5
			atom=Atom(label,(x,y,0.0))
			edge.append(atom)
		unitCell=Atoms()
		for i in range(6):
			newEdge=edge.copy()
			newEdge.rotate('z',i*2*pi/6.0)
			unitCell.extend(newEdge)
		return unitCell
import sys		
if __name__=='__main__':
	#structure().build(int(sys.argv[1]),int(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4])).write()
	pass