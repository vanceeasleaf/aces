from ase import Atoms
from ase import io
def get_cluster(n,i,dis,cluster=[]):
	for j in range(i+1,n):
		if(dis[i,j]<0.001):
			cluster.append(j)
			get_cluster(n,j,dis,cluster)
	return cluster
	
def get_clusters(n,dis):
	#every atom is first
	clusters=[True for i in range(n)]
	for i in range(n):
		if(not clusters[i]):continue
		cluster=get_cluster(n,i,dis,[])
		for j in cluster:
			clusters[j]=False	
			
	return clusters
	
def get_unique_atoms(atoms,mic=True):
	n=len(atoms)
	dis=atoms.get_all_distances(mic=mic)
	clusters=get_clusters(n,dis)
	newatoms=Atoms()
	for i in range(n):
		if clusters[i]:
			newatoms.append(atoms[i])
	
	cell=atoms.get_cell()
	if mic:
		pbc=[1,1,1]
	else:
		pbc=atoms.get_pbc()
	newatoms.set_pbc(pbc)
	newatoms.set_cell(cell)
	if mic:
		newatoms.positions=newatoms.get_positions(wrap=True)
	return newatoms

def atoms_from_dump(filename,elements=None):
	atoms=io.read(filename,format='lammps')	
	if elements:
		s=atoms.numbers
		symbols=[elements[i-1] for i in s ]
		atoms.set_chemical_symbols(symbols)
	return atoms