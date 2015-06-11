from ase import Atoms
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
	
def get_unique_atoms(atoms):
	n=len(atoms)
	dis=atoms.get_all_distances()
	clusters=get_clusters(n,dis)
	newatoms=Atoms()
	for i in range(n):
		if clusters[i]:
			newatoms.append(atoms[i])
	return newatoms