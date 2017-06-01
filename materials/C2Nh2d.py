# encoding : utf8
# C2N hollow 2D structure
from ase import Atoms,Atom
from math import sqrt,pi
from aces import tools
from aces.materials  import Material
import numpy as np
class structure(Material):
	def set_parameters(self):
		self.sideLen=1
		self.cubic=True
	def setup(self):
		pass
	def lmp_structure(self):
		atoms=self.getUnitCell().repeat([self.latx,self.laty,1])
		cell=atoms.cell*self.bond
		#atoms.positions[:,2]+=np.random.uniform(-0.5, 0.5,len(atoms))
		atoms.set_cell(cell,scale_atoms=True)
		atoms.set_pbc([self.xp,self.yp,self.zp])
		atoms.center()
		return atoms
		
		

		
		
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
		#get cell
		dist=(self.sideLen+1)*3.0 #the distance between  2 hole center
		if self.cubic:
			newAtoms=unitCell.copy()
			newAtoms.translate([dist*sqrt(3)/2,dist/2.0,0])
			unitCell.extend(newAtoms)
			unitCell.set_cell([dist*sqrt(3),dist,10.0])
		else:
			cell=np.diag([dist*sqrt(3)/2,dist,10])
			cell[0,1]=dist/2.0
			unitCell.set_cell(cell)
		return unitCell
