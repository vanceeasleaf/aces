from aces.pizza.data import data
from ase import Atoms
import numpy as np
class lammpsdata:
	def __init__(self,atoms,elements=None):
		self.atoms=atoms		
		self.types=self.getTypes()
		if not elements==None:
			self.types=elements
			
	def mergeVec(self,x,y):
		direct=np.cross(x,y)
		if np.allclose(np.linalg.norm(direct),0):
			direct=[0,0,1]
		else:
			direct=direct/np.linalg.norm(direct)
		phi=np.arccos(np.dot(x,y)/np.linalg.norm(x)/np.linalg.norm(y))
		return (direct,phi)
		
	def getTypes(self):
		return list(set(self.atoms.get_chemical_symbols()))
		
	def writedata(self,filename="structure"):
		a=data()
		a.title=self.atoms.get_chemical_formula()
		unit=self.get_rotated_atoms()
		cell=unit.cell
		a.headers['xlo xhi']=[0,cell[0,0]]
		a.headers['ylo yhi']=[0,cell[1,1]]
		a.headers['zlo zhi']=[0,cell[2,2]]
		a.headers['xy xz yz']=[cell[1,0],cell[2,0],cell[2,1]]
		a.headers['atoms']=len(self.atoms)
		a.headers['atom types']=len(self.types)
		atomsdata=[]
		for i,atom in enumerate(unit):
			x=[i+1,self.types.index(atom.symbol)+1]+list(atom.position)
			atomsdata.append(x)
		a.sections['Atoms']=[' '.join(map(str,x))+'\n' for x in atomsdata]
		a.write(filename)
		
	def get_rotated_atoms(self):
		unit=self.atoms.copy()
		direct,phi=self.mergeVec(unit.cell[0],[1,0,0])
		unit.rotate(direct,phi,rotate_cell=True)
		yn=[0,unit.cell[1,1],unit.cell[1,2]]
		direct,phi=self.mergeVec(yn,[0,1,0])
		unit.rotate(direct,phi,rotate_cell=True)
		return unit