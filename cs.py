from aces.runners import Runner
import numpy as np
from aces.tools import *
import h5py
from  numpy import linalg as LA
import os
from ase import io,Atoms
from scipy import optimize
from aces.f import read_forces,writefc2,writefc3,disp2atoms

def shrink(y,a):
	return np.sign(y)*np.maximum(np.abs(y)-a,0.0)
def dot(a,b):
	return np.tensordot(a,b,axes=([1],[0]))
def maxeig(A):
	for j in A.shape[1]:
		B[:,j]=A[:,j]/A.sum(axis=0)[j]
	C=B.sum(axis=1)
	W=C/C.sum()
	lmax=(A.dot(W)/W).average()
	return lmax

class runner:
	def __init__(self,NAH=3,split=True,mu=0.1):
		self.NAH=NAH
		self.split=split
		self.mu=mu

		#self.db=h5py.File('force.hdf5')
	def getForce(self,pos,files):
		print "reading vasprun.xml and POSCAR"
		u=[]
		for file in files:
			dir=os.path.dirname(file)
			atoms=io.read(dir+'/POSCAR')
			u.append(atoms.positions-pos)
		forces=[]
		for file in files:
			forces.append(read_forces(file))
		return np.array(forces),np.array(u)
	def getsatoms(self):
		filename='disp_fc3.yaml'
		if(exists(filename)):
			return disp2atoms(filename)
		filename='disp.yaml'
		if(exists(filename)):
			return disp2atoms(filename)
		filename='3RD.SPOSCAR'
		if(exists(filename)):
			from ase import io
			return io.read(filename,format='vasp')
		filename='SPOSCAR'
		if(exists(filename)):
			from ase import io
			return io.read(filename,format='vasp')

	def getSupercell(self):
		atoms=self.getsatoms()
		from pyspglib import spglib
		s=spglib.get_symmetry(atoms)
		symmetry=[]
		pos=atoms.positions
		print "building symmetry"
		for i,rot in enumerate(s['rotations']):
			trans = s['translations'][i]
			map0=self.getMap(atoms,rot,trans)
			symmetry.append([rot,map0])
		return pos,symmetry
	def getMap(self,atoms,rot,trans):
		v=atoms.copy()
		v.positions=v.positions.dot(rot.T)
		v.translate(trans.dot(v.cell))
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
		return map0
	def getTrainSets(self,u):
		#f=self.db
		#if not 'Fu' in f:
		#	f['Fu']=self.getForce()
		#F,u=f['Fu']

		assert len(u)>0
		self.L=len(u)
		n=self.natom=len(u[0])
		row=0
		rowr=[0]
		for i in range(self.NAH):
			row+=(n*3)**i
			rowr.append(row)
		self.rowr=rowr
		
	def getMatrix(self,F,u):
		self.getTrainSets(u)
		print "getting compressive matrix"
		rowr=self.rowr
		A=np.zeros([self.L,rowr[-1]])
		g=self.mulU
		shape=F.shape
		
		for j in range(self.L):
			for i in range(self.NAH):
				r=range(rowr[i],rowr[i+1])		
				A[j,r]=-g(u[j].flatten(),i)
		F=F.reshape([shape[0],shape[1]*shape[2]])
		return F,A

	def run(self):
		pos,symmetry=self.getSupercell()
		files=shell_exec('find dirs/dir_* -name vasprun.xml |sort -n').split('\n')
		if len(files)>100:	files=files[:100]
		F,u=self.getForce(pos,files)
		if False:
			F0=[]
			u0=[]
			for i,f in enumerate(F):
				for rot,map0 in symmetry:
					newF=f[map0].dot(rot.T)
					newu=u[i][map0].dot(rot.T)
					F0.append(newF)
					u0.append(newu)
			F=np.r_[F,F0]
			u=np.r_[u,u0]
		
		v=LA.norm(u,axis=2)
		u0=v.flatten().max()

		F,A=self.getMatrix(F,u/u0)
		print "start compressive sensing "
		B=cs(mu=self.mu,split=self.split).run(F,A)
		print "rebuilding IFCs "
		phi=self.rebuild(B)
		print "writing IFCs "
		fc2=np.einsum(phi[1],[1,0,3,2])/u0

		writefc2(fc2,'csfc2')
		if self.NAH>=3:
			a=h5py.File('fc3p.hdf5')
			if 'fc3' in a:del a['fc3']
			a['fc3']=phi[2]/u0/u0
			a.close()
			self.fc3()
	def fc3(self):
		print "writing csfc3 "
		a=h5py.File('fc3p.hdf5')
		fc3=np.einsum(a['fc3'],[0,2,1,3,5,4])
		from ase import io
		atoms=io.read('POSCAR')
		satoms=self.getsatoms()
		writefc3(fc3,atoms,satoms,'csfc3')
	def rebuild(self,B):
		#f=self.db
		#if not 'FA' in f:
		#	f['FA']=self.getMatrix()
		#F,A=f['FA']
		
		#print B
		#return B
		n=self.natom
		rowr=self.rowr
		phi=[]
		for i in range(self.NAH):
			r=range(rowr[i],rowr[i+1])
			x=B[r].reshape([n,3]*(i+1))
			idx=np.array([0,i+1])
			rdx=[]
			for j in range(i):
				rdx.extend(idx+(j+1))
			rdx.extend(idx)
			
			x=np.einsum(x,rdx) 
			phi.append(x)

		
		return phi
	
	def mulU(self,x,p):
		if p>0:return np.kron(self.mulU(x,p-1),x)/p
		else: return 1.0
class sym(runner):

	def run(self):
		pos,symmetry=self.getSupercell()
		files=shell_exec('find dirs/dir_* -name vasprun.xml |sort -n').split('\n')
		if len(files)>100:	files=files[:100]
		F,u=self.getForce(pos,files)
		
		v=LA.norm(u,axis=2)
		u0=v.flatten().max()

		F,A=self.getMatrix(F,u/u0,symmetry)
		print "start compressive sensing "
		B=cs(mu=self.mu,split=self.split).run(F,A)
		B=self.restore(B)
		print "rebuilding IFCs "
		phi=self.rebuild(B)
		print "writing IFCs "
		fc2=np.einsum(phi[1],[1,0,3,2])/u0

		writefc2(fc2,'csfc2')
		if self.NAH>=3:
			a=h5py.File('fc3p.hdf5')
			if 'fc3' in a:del a['fc3']
			a['fc3']=phi[2]/u0/u0
			a.close()
			self.fc3()
	def getTrainSets(self,u):
		assert len(u)>0	
		n=self.natom=len(u[0])
		self.L=len(u)
		row=0
		rowr=[0]
		for i in range(self.NAH):
			row+=(n*3)**i
			rowr.append(row)
		self.rowr=rowr
	def getT(self,symmetry):
		n=self.natom
		e=np.eye(3)
		
		for i in range(n):
			for j in range(i):
				x=range(n)
				x[i],x[j]=x[j],x[i]
				symmetry.append([e,x])
		
		return
		Ts=[]
		dim=[[0,0]]
		#build graph		
		for nh in range(self.NAH):
			print "building %d th order graph"%(nh+1)
			G=self.graph(n,symmetry,nh)			
			print "searching for  %d th independent variables"%(nh+1)
			groups=self.connectedGraph(G)
			print "%d th order independent variables: "%(nh+1)
			for g in groups:			
				print self.unpa(n,g[0],nh)
			N=len(groups)*3**(nh+1)
			M=(3*n)**(nh+1)
			T0=np.zeros([M,N])
			dim.append([dim[-1][0]+M,dim[-1][1]+N])
			for i,g in enumerate(groups):
				members=self.getMember(G,nh,g)
				yr=np.arange(i*3**(nh+1),(i+1)*3**(nh+1))
				for id,rot in members:		
					xr=np.arange(id*3**(nh+1),(id+1)*3**(nh+1))
					T0[xr][:,yr]=self.mulR(rot,nh+1)
			Ts.append(T0)
		T=np.zeros(dim[-1])
		for nh in range(self.NAH):
			d1=dim[nh]
			d2=dim[nh+1]
			xr=np.arange(d1[0],d2[0])
			yr=np.arange(d1[1],d2[1])
			T[xr][:,yr]=Ts[nh]
		return T.reshape([3*n,-1,dim[-1][1]])
	def mulR(self,x,p):
		if p>0:return np.kron(self.mulR(x,p-1),x)
		else: return 1.0
	def getMember(self,G,nh,g):		
		c=g[0]
		mem=[[c,np.eye(3)]]
		for r in g[1:]:
			m=self.getMember(G,nh,r[0])
			for x,y in m:
				mem.append([x,y.dot(r[1])])
		return mem
	def connectedGraph(self,G):
		n=len(G)
		occupied=np.zeros(n,dtype=np.bool)
		root=[]
		for i in range(n):
			if occupied[i]:continue
			node=self.tree(G,occupied,i)
			root.append(node)

		return root	
	def tree(self,G,occupied,i):
		root=[i]
		for u in G[i]:
			j=u[0]
			if occupied[j]:continue
			occupied[j]=True
			node=self.tree(G,occupied,j)
			root.append([node,u[1]])
		return root


	def graph(self,natom,symmetry,nh=2):
		N=natom**(nh+1)
		G=[[] for i in range(N)]
		for rot,map1 in symmetry:
			for ii in range(N):
				idx=[map1[i] for i in self.unpa(natom,ii,nh)]
				jj=self.pa(natom,idx,nh)
				G[ii].append([jj,rot])
				G[jj].append([ii,rot.T])
		return G
	def unpa(self,natom,ii,nh=2):
		idx=[]
		for i in range(nh,-1,-1):
			idx.append(ii//natom**i)
			ii=ii%(natom**i)
		return idx
	def pa(self,natom,idx,nh=2):
		s=0
		for i in range(nh,-1,-1):
			s+=idx[nh-i]*natom**i
		return s		
	def getMatrix(self,F,u,symmetry):
		self.getTrainSets(u)
		print "getting compressive matrix"
		T=self.T=self.getT(symmetry)
		self.nIV=T.shape[2]
		rowr=self.rowr
		n=self.natom*3
		A=np.zeros([self.L*n,self.nIV])
		g=self.mulU
		F=F.reshape([self.L*n,1])
		for j in range(self.L):
			a=np.zeros(rowr[-1])
			for i in range(self.NAH):
				r=range(rowr[i],rowr[i+1])	
				a[r]=-g(u[j].flatten(),i)
			for k in range(n):
				A[j+n*k,:]=a.dot(T[k])
		return F,A
	def restore(self,v):
		return np.einsum(T,[1,0,2],v,[2])

class cs:
	def __init__(self,mu=0.7,split=True,lam=0.9):
		self.mu,self.lam=mu,lam
		self.split=split
	def initu(self,f,A):
		dim=list(f.shape)
		dim[0]=A.shape[1]
		#so dim is the shape of u
		return np.ones(dim)
	def testcs(self):
		f=(np.ones(1)*20.0).reshape(1,1)
		A=np.array([7.0,10.0]).reshape(1,2)
		print self.run(f,A)
	def test2(self):
		f=np.array([7.0,8.0])
		A=np.array([[1.0,0],[1.0,0]])
		print self.run(f,A)
	def test3(self):
		f=np.array([7.0,8.0])
		A=np.array([[1.0,0,0],[1.0,0,0]])
		print self.run(f,A)
	def test4(self):
		f=np.array([7.0,8.0]).reshape(1,2)
		A=np.array([[1.0,0]])
		print self.run(f,A)
	def run(self,f,A):

		#normalize 	
		print "normalizing sensing matrix"

		#from scipy.sparse.linalg import eigsh
		"""
		aA=eigsh(A.T.dot(A),k=6)[0].max()#the largest eigenvalue
		f/=np.sqrt(aA)
		#print LA.norm(f)
		A/=np.sqrt(aA)
		#print LA.norm(A)
		"""
		aA=np.double(A.shape[0]**A.max().max())#maxeig(A.T.dot(A))
		f/=np.sqrt(aA)
		A/=np.sqrt(aA)
		"""

		v=np.eye(len(A.T))-A.T.dot(A)
		for i in range(20):
			v=v.dot(v)
			print LA.norm(v)
		"""
		if self.split:return self.split_bregman(f,A)
		else:
			return self.bregman(f,A)

	def split_bregman(self,f,A):
		def cc(u1):
			print "CG error:",LA.norm(u1-self.bb.flatten())/LA.norm(self.bb)
			self.bb=u1
		def g(u,*args):
			A,f,lam,d,mu,b=args
			u=u.reshape(shape)
			return 1.0/2*LA.norm(np.dot(A,u)-f)**2+lam/2.0*LA.norm(d-mu*u-b)**2
		def dg(u,*args):
			A,f,lam,d,mu,b=args
			u=u.reshape(shape)
			return (A.T.dot(A.dot(u)-f)-lam*mu*(d-b-mu*u)).flatten()
		u=self.initu(f,A)
		shape=u.shape
		d=np.zeros_like(u)
		b=np.zeros_like(u)
		deta=0.001
		erru=1.0
		lam=self.lam
		mu=self.mu
		scale=1.0/np.amax(np.abs(f))*1000.0
		print "scale="+str(scale)
		f0=np.zeros_like(f)
		self.bb=np.zeros_like(u)
		#f*=scale
		print 'dimmensions:', A.shape,u.shape
		while erru>deta:			
			#g=lambda u:1.0/2*LA.norm(dot(A,u.reshape(shape))-f)**2+lam/2.0*LA.norm(d-mu*u.reshape(shape)-b)**2
			f1=(f*scale-dot(A,u))+(f0+dot(A,u))/2
			u1=optimize.fmin_cg(g, u,args=(A,f1,lam,d,mu,b),disp=False,fprime=dg,callback=cc,gtol=deta*10).reshape(shape)

			d1=shrink(mu*u1+b,1.0/lam)
			b1=b+mu*u1-d1
			erru=LA.norm(u1-u)/LA.norm(u)
			print 'split bregman iteration error:',erru
			b=b1
			u=u1
			d=d1
			f0=f1
		return u/scale

	def bregman(self,f,A):
		u=self.initu(f,A)
		f0=np.zeros_like(f)
		deta=0.0001
		erru=1
		scale=1000.0
		while erru>deta:			
			f1=f*scale+f0-dot(A,u)
			u1=self.FCP(f1,A,u)
			
			erru=LA.norm(u1-u)/LA.norm(u)
			print 'bregman iteration error:',erru
			u=u1
			f0=f1
		return u/scale

	def FCP(self,f,A,u=None):
		if u is None:
			u=self.initu(f,A)

		m,n=A.shape
		if len(f.shape)>1:
			n*=list(f.shape)[1]
		ta=1.99999#min(1.999,max(1.1,-1.665*np.float(m)/n+2.665))
		
		mu=self.mu
		deta=0.01
		errg=1
		erru=1
		while  erru>deta :#or errg>deta:
			p=np.dot(A,u)-f
			
			g=np.dot(A.T,p)
			
			u1=shrink(u-ta*g,mu*ta)
			errg=1.0/mu*LA.norm(g,np.inf)-1
			erru=LA.norm(u1-u)/LA.norm(u)
			print 'FCP iteration :',erru
			u=u1
			
		return u
class energy(runner):
	def run(self):
		pos,symmetry=self.getSupercell()
		files=shell_exec('find dirs/dir_* -name vasprun.xml |sort -n').split('\n')	
		F,u=self.getForce(pos,files)
		V=(F*u).sum(axis=(1,2))
		v=LA.norm(u,axis=2)
		u0=v.flatten().max()
		F,A=self.getMatrix(V,u/u0)
		print "start compressive sensing "
		B=cs(mu=self.mu,split=self.split).run(F,A)
		print "rebuilding IFCs "
		phi=self.rebuild(B)
		print "writing IFCs "
		fc2=np.einsum(phi[1],[1,0,3,2])/u0

		writefc2(fc2,'csfc2')
		if self.NAH>=3:
			a=h5py.File('fc3p.hdf5')
			if 'fc3' in a:del a['fc3']
			a['fc3']=phi[2]/u0/u0
			a.close()
			self.fc3()
	def getTrainSets(self,u):
		#f=self.db
		#if not 'Fu' in f:
		#	f['Fu']=self.getForce()
		#F,u=f['Fu']

		assert len(u)>0
		self.L=len(u)
		n=self.natom=len(u[0])
		row=0
		rowr=[0]
		for i in range(self.NAH):
			row+=(n*3)**(i+1)
			rowr.append(row)
		self.rowr=rowr
		
	def getMatrix(self,F,u):
		self.getTrainSets(u)
		print "getting compressive matrix"
		rowr=self.rowr
		A=np.zeros([self.L,rowr[-1]])
		g=self.mulU
		shape=F.shape
		
		for j in range(self.L):
			for i in range(self.NAH):
				r=range(rowr[i],rowr[i+1])		
				A[j,r]=-g(u[j].flatten(),i+1)
		return F,A