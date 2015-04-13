# encoding : utf8
# C2N hollow 2D structure
from ase import Atoms,Atom
from math import sqrt,pi
from ase.io.vasp import write_vasp
from ase.UnitCell.unitcell import UnitCell
#import sys
#sys.path.append("/home1/xggong/zhouy/tcscripts/pizza/")
#from data import data
class C2Nh2D:
	def __init__(self,latx,laty,bond,sideLen):
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
	
	def write(self):
		self.atoms.write("CN.xyz")
		write_vasp("POSCAR",self.atoms,sort="True",direct=True,vasp5=True)
		poscar = open("POSCAR")
		unit_cell = UnitCell(poscar)
		unit_cell.num_atom_types=2
		lammps=open("structure","w")
		lammps.write(unit_cell.output_lammps())
		lammps.close()
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
	C2Nh2D(int(sys.argv[1]),int(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4])).write()
	