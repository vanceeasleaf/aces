from aces.pizza.data import data
from ase import Atoms,Atom
import numpy as np
class lammpsdata:
	def __init__(self,atoms=None,elements=None):
		if atoms is None:
			atoms=Atoms()
		self.atoms=atoms		
		self.types=self.getTypes()
		if not elements is None:
			self.types=elements
	def set_src(self,filename):
		a=data(filename)	
		a.map(1,'id',2,'type',3,'x',4,'y',5,'z')
		ats=a.viz(0)[2]	
		atoms=Atoms()
		for at in ats:
			id,type,x,y,z=at
			if len(self.types)==0:
				atoms.append(Atom(type,position=[x,y,z]))
			else:
				atoms.append(Atom(self.types[type-1],position=[x,y,z]))
		cell=np.zeros([3,3])
		xlo,xhi=a.headers["xlo xhi"]
		ylo,yhi=a.headers["ylo yhi"]
		zlo,zhi=a.headers["zlo zhi"]
		
		cell[0,0]=xhi-xlo
		cell[1,1]=yhi-ylo
		cell[2,2]=zhi-zlo
		if a.headers.has_key("xy xz yz"):
			xy,xz,yz=a.headers["xy xz yz"]
			cell[1,0],cell[2,0],cell[2,1]=xy,xz,yz
		atoms.set_cell(cell)
		self.atoms=atoms
		return atoms
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
		unit,rot=self.get_rotated_atoms()
		cell=unit.cell
		a.headers['xlo xhi']=[0,cell[0,0]]
		a.headers['ylo yhi']=[0,cell[1,1]]
		a.headers['zlo zhi']=[0,cell[2,2]]
		#if(cell[1,0]>cell[0,0]*0.5):cell[1,0]=cell[0,0]/abs(cell[0,0])*int((abs(cell[0,0])*0.5)*1000000)/1000000.0
		v=1.0
		if(cell[1,0]<0):v=-1.0
		cell[1,0]=v*int((abs(cell[1,0]))*10000)/10000.0
		#print [cell[1,0]/cell[0,0],cell[2,0]/cell[0,0],cell[2,1]/cell[1,1]];
		if not np.allclose([0,0,0],[cell[1,0]/cell[0,0],cell[2,0]/cell[0,0],cell[2,1]/cell[1,1]],atol=0.01):
			a.headers['xy xz yz']=[cell[1,0],cell[2,0],cell[2,1]]
		a.headers['atoms']=len(self.atoms)
		a.headers['atom types']=len(self.types)
		atomsdata=[]
		for i,atom in enumerate(unit):
			x=[i+1,self.types.index(atom.symbol)+1]+list(atom.position)
			atomsdata.append(x)
		a.sections['Atoms']=[' '.join(map(str,x))+'\n' for x in atomsdata]
		a.write(filename)
		return rot
	def get_rotated_atoms(self):
		unit=self.atoms.copy()
		direct,phi=self.mergeVec(unit.cell[0],[1,0,0])
		unit.rotate(direct,phi,rotate_cell=True)
		yn=[0,unit.cell[1,1],unit.cell[1,2]]
		direct1,phi1=self.mergeVec(yn,[0,1,0])
		unit.rotate(direct1,phi1,rotate_cell=True)
		return unit,(direct,phi,direct1,phi1)