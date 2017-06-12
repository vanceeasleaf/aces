from ase import Atoms
from ase import io
import numpy as np
# rotate atoms to swap the x z axis for fix=1 and so on, keep the axis right hand
def swap(atoms,fix=1):
	direct=[1,1,1]
	direct[fix]=0
	atoms.rotate(direct,np.pi,rotate_cell=True)
	order=[[0,2,1],[2,1,0],[1,0,2]][fix]
	cell=atoms.cell[order]
	cell[fix]*=-1
	atoms.set_cell(cell)	
	
def center_box(atoms):
	offset=np.sum(atoms.cell,axis=0)/2
	atoms.translate(-offset)
def extent(atoms):
	return atoms.positions.max(axis=0)-atoms.positions.min(axis=0)

def getMassFromLabel(labels):
	from ase.data import atomic_masses,atomic_numbers
	nums=[atomic_numbers[label] for label in labels]
	masses=[atomic_masses[num] for num in nums]
	return masses
def get_cluster(n,i,dis,cluster=[]):
	for j in range(i+1,n):
		if(dis[i,j]<0.01):
			cluster.append(j)
			get_cluster(n,j,dis,cluster)
	return cluster
	
def get_clusters(n,dis):
	"""find the head of link-table made of atoms at the same positions
	
	[description]
	
	Arguments:
		n {[type]} -- [number of atoms]
		dis {[array[n][n]} -- [distance between atoms]
	
	Returns:
		[list] -- [the boolean value wheather that atom is at the same position of some other atom]
	"""
	#every atom is head
	clusters=[True]*n
	for i in range(n):
		# if the position of atom is occupied
		if(not clusters[i]):continue
		cluster=get_cluster(n,i,dis,[])
		for j in cluster:
			clusters[j]=False	
			
	return clusters
def wrap(atoms):
	atoms.positions=atoms.get_positions(wrap=True)

def get_unique_atoms(atoms,mic=True):
	n=len(atoms)
	u=atoms.copy()
	if(mic):
		u.set_pbc([True]*3)
	dis=u.get_all_distances(mic=mic)
	clusters=get_clusters(n,dis)
	newatoms=Atoms()
	for i in range(n):
		if clusters[i]:
			newatoms.append(atoms[i])
	
	cell=atoms.get_cell()
	pbc=atoms.get_pbc()
	newatoms.set_pbc(pbc)
	newatoms.set_cell(cell)
	#newatoms.center()
	return newatoms

def atoms_from_dump(filename,elements=None,index=-1):
	atoms=io.read(filename,format='lammps',index=index)	
	if elements:
		s=atoms.numbers
		symbols=[elements[i-1] for i in s ]
		atoms.set_chemical_symbols(symbols)
	return atoms