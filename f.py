import numpy as np
from aces.tools import *
def writefc2(fc2,filename='FORCE_CONSTANTS_2ND'):
	natom=len(fc2)
	s="%d\n"%natom
	for i in range(natom):
		for j in range(natom):
			s+="%d\t%d\n"%(i+1,j+1)
			s+=matrixFormat(fc2[i,j])
	write(s,filename)
def writefc3(fc3,atoms,satoms,filename='FORCE_CONSTANTS_3RD'):
	psm=premitiveSuperMapper(atoms,satoms)
	natom=len(fc3)
	s=""
	n=0
	scale=np.amax(np.abs(fc3))
	for i in range(natom):
		for j in range(natom):
			for k in range(natom):
				if np.allclose(np.zeros([3,3,3]),fc3[i,j,k]/scale,atol=1e-04):continue
				n+=1
				s+="\n%d\n"%n
				s+=psm.maps2p(i,j,k)
				s+=matrix3Format(fc3[i,j,k])
	s="%d\n"%n+s			
	write(s,filename)
class premitiveSuperMapper:
	def __init__(self,atoms,satoms):
		self.atoms=atoms
		self.satoms=satoms
	def getS2p(self):
		atoms=self.atoms
		cellp=atoms.cell
		posp=atoms.get_scaled_positions()
		satoms=self.satoms
		cell=satoms.cell
		pos=satoms.get_scaled_positions()

		#which cell is i j in
		vpos=pos*np.linalg.norm(cell,axis=1)/np.linalg.norm(cellp,axis=1)
		v=np.floor(vpos)
		
		celloffset=v.dot(cellp)
		vpos-=v
		s2p=-np.ones(len(vpos),dtype='int')
		for pd,p in enumerate(vpos):
			for id,a in enumerate(posp):
				if np.allclose(p,a):
					s2p[pd]=id+1
					break

		return s2p,celloffset
	def maps2p(self,i,j,k):
		if not hasattr(self, 's2p'):
			self.s2p,self.celloffset=self.getS2p()
		v=self.celloffset

		s=toString(v[j]-v[i])+'\n'
		s+=toString(v[k]-v[i])+'\n'
		s+='%d %d %d\n'%(self.s2p[i],self.s2p[j],self.s2p[k])
		return s

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
def matrixFormat(mat):
	n,m=mat.shape
	s=""
	for k in range(n):
		for p in range(m):
			s+="\t%f"%mat[k,p]
		s+="\n"	
	return s
def matrix3Format(mat):
	n,m,q=mat.shape
	s=""
	for k in range(n):
		for p in range(m):
			for t in range(q):
				s+=" %d  %d  %d       %f\n"%(k+1,p+1,t+1,mat[k,p,t])
	return s
def rotationMatrix(axis,theta):
	from numpy import cross,eye,dot
	from scipy.linalg import expm3,norm
	return expm3(cross(eye(3),axis/norm(axis)*theta))
def RotateVector(vec,axis,theta):
	#debug(rotationMatrix(axis,theta))
	return np.dot(rotationMatrix(axis,theta),vec)
def toString(m,sep=' '):
	return sep.join(map(str,m))
def refinefc3():
	f=open('FORCE_CONSTANTS_3RD')
	g=open('fc3new','w')
	nblock=int(f.next().split()[0])
	fc3=[]
	r1t=[]
	r2t=[]
	idx=[]
	#print >>g,nblock
	for i in range(nblock):
		f.next()
		f.next()
		r1=np.array(map(float,f.next().split()))
		r1t.append(r1)
		r2=np.array(map(float,f.next().split()))
		r2t.append(r2)
		idx.append(map(int,f.next().split()))
		fc=np.zeros([3,3,3])
		for i in range(3):
			for j in range(3):
				for k in range(3):
					fc[i,j,k]=float(f.next().split()[3])
		fc3.append(fc)
	scale=np.amax(np.abs(fc3))
	u=[]
	for i in range(nblock):
		if np.allclose(np.zeros([3,3,3]),fc3[i]/scale,atol=1e-01*.3) :
			u.append(False)
		else:
			u.append(True)
	u=np.array(u)
	fc3=np.array(fc3)[u]
	r1t=np.array(r1t)[u]
	r2t=np.array(r2t)[u]
	idx=np.array(idx)[u]
	v=np.unique(idx[:,0])
	filters=np.array([idx[:,0]==c for c in v])
	
	for i,f in enumerate(filters):
		a=np.arange(len(idx))[f][0]
		fc3[a]-= fc3[f].sum(axis=0)	
	v=np.unique(idx[:,1])
	filters=np.array([idx[:,1]==c for c in v])
	
	for i,f in enumerate(filters):
		a=np.arange(len(idx))[f][0]
		fc3[a]-= fc3[f].sum(axis=0)	
	nblock=len(fc3)
	print >>g,nblock
	for i in range(nblock):
		print >>g,i+1
		print >>g,toString(np.around(r1t[i],3))
		print >>g,toString(np.around(r2t[i],3))
		print >>g,toString(idx[i])
		print >>g,matrix3Format(fc3[i])
def rotatefc3(t):
	f=open('FORCE_CONSTANTS_3RD')
	g=open('fc3new','w')
	M=rotationMatrix([0,0,1],t*np.pi/180.0)
	nblock=int(f.next().split()[0])
	print >>g,nblock
	for i in range(nblock):
		f.next()
		print >>g,f.next()
		r1=np.array(map(float,f.next().split()))
		r1t=M.dot(r1)
		print >>g,toString(r1t)
		r1=np.array(map(float,f.next().split()))
		r1t=M.dot(r1)
		print >>g,toString(r1t)
		print >>g,f.next()
		fc=np.zeros([3,3,3])
		for i in range(3):
			for j in range(3):
				for k in range(3):
					fc[i,j,k]=float(f.next().split()[3])

		fc=np.einsum('ij,jkl',M,fc)
		fc=np.einsum('ij,kjl',M,fc)
		fc=np.einsum('ij,klj',M,fc)
		print >>g,matrix3Format(fc)

try:
    from lxml import etree as ElementTree
    xmllib="lxml.etree"
except ImportError:
    try:
        import xml.etree.cElementTree as ElementTree
        xmllib="cElementTree"
    except ImportError:
        import xml.etree.ElementTree as ElementTree
        xmllib="ElementTree"
def writevasp(atoms,file='POSCAR'):
	f=open(file,'w')
	s=np.array(atoms.get_chemical_symbols())
	ss=atoms.get_scaled_positions()
	print >>f,'ACES POSCAR'
	print >>f,'1.0'
	for x in atoms.cell:
		print >>f,toString(x)
	#ele=np.unique(s)
	ele=[]
	for a in s:
		if a in ele:
			continue
		ele.append(a)
	print >>f,toString(ele)
	a=[]
	p=np.arange(len(s))
	for e in ele:
		a.append(p[s==e])
	ns=[len(x) for x in a]
	print >>f,toString(ns)
	print >>f,'Direct'
	v=[]
	for x in a:
		for u in x:
			v.append(u)
			print >>f,toString(ss[u])
	f.close()
	x=np.array(v,dtype=np.int).argsort()
	np.savetxt('POSCARswap',x)
def mapatoms(newatoms,oldatoms):
	#the structure is the same but order changed , return the index ,so that newatoms[index]=oldatoms
	atoms=newatoms
	v=oldatoms
	import itertools
	from scipy.spatial.distance import cdist
	posi=atoms.positions
	d2s=np.empty((27,len(v),len(v)))
	for j,(ja,jb,jc) in enumerate(itertools.product(xrange(-1,2),
                                                xrange(-1,2),
                                                xrange(-1,2))):
		posj=v.positions+np.dot([ja,jb,jc],v.cell)
		d2s[j,:,:]=cdist(posi,posj,"sqeuclidean")
	d2min=d2s.min(axis=0)
	map0=np.argmin(d2min,axis=1)
	return atoms[map0].positions,map0
	#return ktmatch1(oldatoms.positions,newatoms.positions)
def ktmatch1(a,b):
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
def read_forces(filename):
    """
    Read a set of forces on atoms from filename, presumably in
    vasprun.xml format.
    """
    calculation=ElementTree.parse(filename
                                  ).getroot().find("calculation")
    for a in calculation.findall("varray"):
        if a.attrib["name"]=="forces":
            break
    nruter=[]
    for i in a.getchildren():
        nruter.append([float(j) for j in i.text.split()])
    nruter=np.array(nruter,dtype=np.double)
    return nruter
def readfc2(filename='FORCE_CONSTANTS'):
		f=open(filename)
		line=f.next()
		natom=int(line)
		fc=np.zeros([natom,natom,3,3])
		for i in range(natom):
			for j in range(natom):
				f.next()
				for k in range(3):
					fc[i,j,k]=map(float,f.next().split())
		return fc