import numpy as np
from aces.tools import *
from numpy.linalg  import norm
import time
from functools import wraps
np.set_printoptions(precision=3,suppress=True)
def fn_timer(function):
	@wraps(function)
	def function_timer(*args, **kwargs):
		t0 = time.time()
		result = function(*args, **kwargs)
		t1 = time.time()
		print ("Total time running %s: %s seconds" %
			(function.func_name, str(t1-t0))
			)
		return result
	return function_timer


def disp2atoms(disp='disp.yaml'):
	from aces.tools import parseyaml
	data = parseyaml(disp)
	cell=np.array(data['lattice'])
	satoms=data['atoms']
	symbols=''.join([o['symbol'] for o in satoms])
	pos=np.array([o['position'] for o in satoms])
	from ase import Atoms
	atoms=Atoms(symbols,scaled_positions=pos,cell=cell)
	return atoms

def rotationMatrix(axis,theta):
	from numpy import cross,eye,dot
	from scipy.linalg import expm3,norm
	return expm3(cross(eye(3),axis/norm(axis)*theta))
def RotateVector(vec,axis,theta):
	#debug(rotationMatrix(axis,theta))
	return np.dot(rotationMatrix(axis,theta),vec)


def others_match(unit,supercell):
	"""find a given unit from a super cell by moving and rotating the unit cell
	
	return a mapping map0 so that fc_s=fc[map0][:,map0].
	fc is of unit,fc_s is of supercell
	atom i in supercell is translated from atom j in unit
	map[i]=j
	fc_s[a,b]=fc[map0][:,map0][a,b]=fc[map0[a]][map0[b]]

	don't forget the atom type must match
	supercell.cell=np.einsum('im,ml->il',rot,unit.cell)=np.dot(rot,unit.cell)

	forces rotate as unit.cell
	fc_s1=np.einsum('im,ml->il',rot,fc_s)
	fc_s1=np.einsum('im,km->ik',rot,fc_s1)
	=>fc_s1=rot*fc_s*rot.T
	supercell.cell[0]==i*np.dot(rot,unit.cell)[0]
	supercell.cell[1]==j*np.dot(rot,unit.cell)[1]
	supercell.cell[2]==k*np.dot(rot,unit.cell)[2]

	Arguments:
		unit {[type]} -- [description]
		supercell {[type]} -- [description]
		return map0,rot
	"""
	pass
	unit=unit.copy()
	supercell=supercell.copy()
	t=find_tripple(supercell,unit)
	#t is the found index of first 3 atoms of unit in supercell
	assert t is not None
	#We are going to find the transform to move unit to the target 3 atoms
	tatoms=supercell[t]
	unit,rot=find_transform(unit,tatoms)
	map0=mapatoms(supercell,unit)
	direct,phi,direct1,phi1=rot
	M1=rotationMatrix(direct,phi)
	M2=rotationMatrix(direct1,phi1)
	#r'=M2.M1.r
	rot=M2.dot(M1)
	print rot
	return map0,rot

def find_transform(unit,tatoms):
	"""find the transform to move unit to tatoms
	
	tranlate:unit[0]->tatoms[0]
	then rotate unit[1]->tatoms[1]
	then rotate unit[2]->tatoms[2]
	Arguments:
		unit {Atoms[3]} -- [description]
		tatoms {Atoms[3]} -- [description]
	"""
	unit=unit.copy()
	tranlate=tatoms.positions[0]-unit.positions[0]
	unit.translate(tranlate)
	direct,phi=merge_vector(unit.positions[1],tatoms.positions[1])
	unit.rotate(direct,phi,rotate_cell=True)
	direct1,phi1=merge_vector(unit.positions[2],tatoms.positions[2])
	return unit,(direct,phi,direct1,phi1)

def merge_vector(x,y):
	"""find the direction and rotation angle to merge vector x to y
	
	[description]
	
	Arguments:
		x {array(3)} -- the operated vector
		y {array(3)} -- the target vector
	"""
	# the rotation dirction is vertical both to x and y
	 
	
	# getting unit direction vector
	direct=np.cross(x,y)
	if np.allclose(np.linalg.norm(direct),0):
		direct=[0,0,1]
	else:

		direct=direct/np.linalg.norm(direct)

	# find the angle between x and y	
	phi=np.arccos(np.dot(x,y)/np.linalg.norm(x)/np.linalg.norm(y))
	return (direct,phi)	

def find_tripple(supercell,unit):
	pos=supercell.positions
	posu=unit.positions
	sys=supercell.get_chemical_symbols()
	sysu=unit.get_chemical_symbols()
	N=len(supercell)
	err=0.01
	for i in xrange(N):
		if sys[i]!=sysu[0]:
			continue
		ox=pos[i]-posu[0]
		for j in xrange(N):
			if norm(pos[j]-posu[1]-ox)>err or sys[j]!=sysu[1]:
				continue
			for k in xrange(N):
				if norm(pos[k]-posu[2]-ox)>err or sys[k]!=sysu[2]:
					continue
				return [i,j,k]
	return None
def selfmapping(newatoms,oldatoms):
	return mapatoms(newatoms,oldatoms)

def mapforce(fc,map0):
	"""when the order of atoms change,please give the corresponding force constants
	
	force of atom i in oldatoms(supercell)->cell rotate ->force of atom i in newatoms(supercell)->wrap to original(supercell) ->force of atom j in original
	we can get mapping as map0, so what's the new order of fc
	
	once rotated we elementwisely rotate fc[i,i1]->fc_rotated[i,i1]

	after wrapping the order must be corrected to use the order of original, in order to compare the rotated fc and the original fc directly.
	
	consider a 1d chain, whose period=3

	|-O-k-O-k-O-k|-O-k-O-k-O-k|-O-k-O-k-O-k|-O-k-O-k-O-k
     the next-neighbor interaction is f
	the fc is    
	-2k-2f     k+f       k+f
	k+f      -2k-2f    k+f
	k+f       k+f      -2k-2f

	use anathor period =4
	|-O-k-O-k-O-k-O-k|-O-k-O-k-O-k-O-k|-O-k-O-k-O-k-O-k|-O-k-O-k-O-k-O-k
	the fc is    
	-2k-2f     k     2f    k
	k     -2k-2f      k     f
	2f       k       -2k-2f   k
	k       2f       k     -2k-2f
	
	    
	to generate fc(4) from fc(3)
	firstly directly repeat the matrix,fc[3,i]=fc[0,i] and fc[i,3]=fc[i,0]
	-2k-2f     k+f    k+f    -2k-2f
	k+f      -2k-2f    k+f    k+f  
	k+f        k+f   -2k-2f   k+f  
	-2k-2f     k+f    k+f    -2k-2f 
	however f and k are undistinguashable just from the fc,assume the cell is large enough 
	Arguments:
		fc {np.array_2d[N,N]} -- [description]
		map0 {np.array_1d[M]->0~N} -- [description]
		return {np.array_2d[M,M]}
	"""
	return fc[map0][:,map0]
def mapatoms(newatoms,oldatoms):

	"""
	atom i in oldatoms-> cell rotate->atom i in newatoms->wrap to unitcell->atom j in unitcell
	so newatoms[i]=oldatoms[j];set(newatoms)===set(oldatoms)*n

	newatoms are wrapped supercell with len(newatoms)=n*len(oldatoms) but multi newatoms->one oldatoms at the same position

	cell rotate could be extended to translation or swapping.
	the structure is the same but order changed , return the mapping ,so that newatoms==oldatoms[map0]
	map0[i]=j=>newatoms[i]==oldatoms[map0[i]]
	len(map0)=len(newatoms)
	oldatoms[map0]=[oldatoms[map0[j]] for j,a in enumerate(map0)]=[newatoms[j] for j,a in enumerate(map0)]=newatoms
	we also have oldatoms[map0[j]]=oldatoms[map0][j]

	rotation dont't change order ,but wrap change
	Arguments:
		newatoms {Atoms[n*N]} -- [description]
		oldatoms {Atoms[N]} -- [description]
	"""		
	atoms=newatoms.copy()
	atoms.cell=oldatoms.cell
	atoms.set_positions(atoms.get_positions(wrap=True))
	v=oldatoms
	import itertools
	from scipy.spatial.distance import cdist
	posi=atoms.positions
	d2s=np.empty((27,len(atoms),len(oldatoms))) #d2s[j,ii,jj] means the distance of atom ii and atom jj with ii in unitcell but jj in translated[ja,jb,jc] unitcell
	for j,(ja,jb,jc) in enumerate(itertools.product(xrange(-1,2),
                                                xrange(-1,2),
                                                xrange(-1,2))):
		posj=v.positions+np.dot([ja,jb,jc],v.cell)
		d2s[j,:,:]=cdist(posi,posj,"sqeuclidean")
	d2min=d2s.min(axis=0) #d2min[ii,jj] means the nearest distance of atom ii and atom jj consider for all the translation
	map0=np.argmin(d2min,axis=1) 
	#print atoms.get_scaled_positions()-oldatoms.get_scaled_positions()[map0]
	print "mapatoms-distances:",np.array([d2min[i,map0[i]] for i in range(len(newatoms))])
	return map0
	#return ktmatch1(oldatoms.positions,newatoms.positions)
""" special KT match
	a=[1,3.1,5.1,5.1]
	b=[3,5,7]
	return order=[null,0,1,null] with tolerance=0.5
"""

def ktmatch1(a,b):
	""" KT match
	@a=[1,3,5,7,6,9,2,8]
	@b=[2,5,3,6,9,8,1,7]
	return @order with order[0]=1 (find 1 at i=6),order[1]=2 (find 3 at i=2)
	order=[b.index(i) for i in a]
	"""
	a=np.array(a)
	b=np.array(b)
	n=len(a)
	order=np.zeros(n,dtype=np.int)
	occupied=np.zeros(n,dtype=np.int)
	for i in range(n):
		s=1000000.0
		jmin=0
		for j in range(n):
			if occupied[j]:continue
			if s>dis(a[i],b[j]):
				s=dis(a[i],b[j])
				jmin=j
		occupied[jmin]=True
		order[i]=jmin
	return b[order],order
def dis(x,y):
	return np.linalg.norm(x-y)		
def ktmatch(a,b):
	a=np.array(a)
	b=np.array(b)
	n=len(a)
	order=range(n)
	if n==1:return a,order
	u,o=ktmatch(a[:-1],b[:-1])	
	
	p=np.zeros_like(a)
	p[:-1]=u
	p[-1]=b[-1]
	neworder=np.r_[o,order[-1]]
	d=dis(a,p)
	for i in range(n-1):
		p[i],p[-1]=p[-1],p[i]
		if d<dis(a,p):
			p[i],p[-1]=p[-1],p[i]
			continue
		neworder[i],neworder[-1]=neworder[-1],neworder[i]
		d=dis(a,p)
	return p,neworder	

