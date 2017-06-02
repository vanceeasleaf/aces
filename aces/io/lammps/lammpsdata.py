from aces.libs.pizza.data import data
from ase import Atoms,Atom
import numpy as np
from aces.tools import debug,toString
from aces.f import merge_vector
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
		bond=False
		if a.headers.has_key("bonds"):
			bond=True
		if bond:
			a.map(1,'id',2,'type',3,'mtype',4,'x',5,'y',6,'z')
		else:
			a.map(1,'id',2,'type',3,'x',4,'y',5,'z')
		
		ats=a.viz(0)[2]	
		#print ats
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

	def getTypes(self):
		return list(set(self.atoms.get_chemical_symbols()))
		
	def writedata(self,filename="structure",creatbonds=-1.0):
		a=data()
		a.title=self.atoms.get_chemical_formula()
		unit,rot=self.get_rotated_atoms()
		cell=unit.cell
		a.headers['xlo xhi']=[0,cell[0,0]]
		a.headers['ylo yhi']=[0,cell[1,1]]
		a.headers['zlo zhi']=[0,cell[2,2]]
		if(cell[1,0]>cell[0,0]*0.5):cell[1,0]=cell[0,0]*0.5
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
			if creatbonds>0.0:
					x=[i+1,1,self.types.index(atom.symbol)+1]+list(atom.position)
			atomsdata.append(x)
		a.sections['Atoms']=[toString(x)+'\n' for x in atomsdata]
		if creatbonds>0.0:
			dis=unit.get_all_distances(mic=True)
			bonds=[]

			n=0
			for i in xrange(len(unit)):
				for j in xrange(i):
					if dis[i,j]<creatbonds:
						n+=1
						bonds.append([n,1,i+1,j+1])
			a.sections['Bonds']=[toString(x)+'\n' for x in bonds]	
			a.headers['bonds']=len(bonds)
			a.headers['bond types']=1

		a.write(filename)
		return rot
	def get_rotated_atoms(self):
		unit=self.atoms.copy()
		#debug(unit.cell)
		direct,phi=merge_vector(unit.cell[0],[1,0,0])
		unit.rotate(direct,phi,rotate_cell=True)
		#debug(unit.cell)
		yn=[0,unit.cell[1,1],unit.cell[1,2]]
		direct1,phi1=merge_vector(yn,[0,1,0])
		unit.rotate(direct1,phi1,rotate_cell=True)
		#debug(unit.cell)
		return unit,(direct,phi,direct1,phi1)