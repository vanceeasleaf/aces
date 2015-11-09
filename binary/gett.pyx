import numpy as np

def gettnh(n,int[:,:] maps,double[:,:,:] rots,nh):

	print "building %d th order graph"%(nh+1)
	G=graph(n,maps,nh)			
	print "searching for  %d th independent variables"%(nh+1)
	groups=connectedGraph(G)
	print "%d th order independent variables: "%(nh+1)
	#for g in groups:			
	#	print unpa(n,g[0],nh)
	N=len(groups)*3**(nh+1)
	M=(3*n)**(nh+1)
	T0=np.zeros([M,N])
	for i,g in enumerate(groups):
		members=getMember(G,nh,g)
		yr=np.arange(i*3**(nh+1),(i+1)*3**(nh+1))
		for id,rot in members:		
			xr=np.arange(id*3**(nh+1),(id+1)*3**(nh+1))
			T0[xr][:,yr]=mulR(rot,nh+1)
	return T0,N
def mulR(x,p):
	if p>0:return np.kron(mulR(x,p-1),x)
	else: return 1.0

def getMember(G,nh,g):		
	c=g[0]
	mem=[[c,np.eye(3)]]
	for r in g[1:]:
		m=getMember(G,nh,r[0])
		for x,y in m:
			mem.append([x,y.dot(r[1])])
	return mem
cdef connectedGraph(int[:,:] G):
	cdef int n=len(G)
	cdef int[:] occupied=np.zeros(n,dtype=np.bool)
	root=[]
	for i in range(n):
		if occupied[i]:continue
		node=tree(G,occupied,i)
		root.append(node)

	return root	
def tree(G,occupied,i):
	root=[i]
	for u in G[i]:
		j=u[0]
		if occupied[j]:continue
		occupied[j]=True
		node=tree(G,occupied,j)
		root.append([node,u[1]])
	return root
def symMatrix():
	M=basicMatrix()
	N=paer0.N
	G=tree()
	headers=np.unique(G[:,1])
	h=[] 
	for j in headers:
		h.append([])
	for i in xrange(N):
		u=G[i,1]
		h[u].append(i)
	P=[[[]for j in headers] for i in range(paer0.natom)]
	for i in range(natom):
		for j in xrange(N):
			x=G[j,1]
			P[i,x].append([M[i,x,0],M[i,x,1],G[j,2]])
	

cdef int[:,:,:] basicMatrix(paer paer0):
	cdef int natom=paer0.natom
	cdef int nh=paer0.nh
	cdef int N=paer0.N
	cdef int m
	cdef int[:,:,:] M=-np.one([natom,N,nh],dtype=np.intc)
	cdef int[:] idx=np.zeros(nh+1,dtype=np.intc)
	for i in xrange(natom):
		for j in xrange(natom):
			for k in xrange(j+1):
				idx[0]=i
				idx[1]=j
				idx[2]=k
				m=paer0.pa(idx)
				M[i,m,0]=j
				M[i,m,1]=k
	return M
cdef int[:,:] tree(int natom,int[:,:] maps,int nh):
	paer0=paer(natom,nh)
	cdef int N=paer0.N
	print "N=%d coefficients,building tree"%N
	print "Number of symmetry operations=%d"%len(maps)
	#3 collunms are parent  , header,rotation respectively
	cdef int[:,:] G=np.zeros([N,3],dtype=np.intc)
	cdef int[:,:] child=-np.ones(N,dtype=np.intc)
	G[:,0]=np.arange(N)
	G[:,1]=np.arange(N)	
	for i,map1 in enumerate(maps):
		for ii in xrange(N):	
			#the element has been some one's children
			if not G[ii,0]==ii:continue	
			jj=paer0.map(ii,map1)
			#these two elements has been the same group and a circle arises
			if G[ii,1]==G[jj,1]:continue
			#set parent
			G[ii,0]=jj
			#set child
			child[jj]=ii
			#set header
			h=ii
			while not h==-1:
				G[h,1]=G[jj,1]
				h=child[h]
			p=jj
			while not child[p]==-1:
				p=child[p]
			child[p]=ii
			#set rotation
			G[ii,2]=i
	return G
cdef extern from "math.h":
	extern long pow(int i,int j) 
cdef class paer:
	cdef int[:] natompow,Ns
	cdef int natom,nh,N
	cdef int[:,:] idxs,idx1
	def __cinit__(self,int natom,int nh):
		self.natom=natom
		self.nh=nh
		
		cdef int[:] natompow=np.zeros(nh+1,dtype=np.intc)
		for i in xrange(nh,-1,-1):
			natompow[i]=pow(natom,i)
		self.natompow=natompow
		self.idx=np.zeros(nh+1,dtype=np.intc)
		self.idx1=np.zeros(nh+1,dtype=np.intc)
		self.idxs=np.zeros([natom**(nh+1),nh+1],dtype=np.intc)
		self.Ns=np.zeros(natom**(nh+1),dtype=np.intc)
		cdef int n=0,m
		if nh==0:
			for i in xrange(natom):
				self.idx[0]=self.idxs[n,0]=i
				m=self.pa(self.idx)				
				self.Ns[m]=n
				n+=1
		if nh==1:
			for i in xrange(natom):
				for j in xrange(i+1):
					self.idx[0]=self.idxs[n,0]=i
					self.idx[1]=self.idxs[n,1]=j
					m=self.pa(self.idx)				
					self.Ns[m]=n
					n+=1			
		if nh==2:			
			for i in xrange(natom):
				for j in xrange(i+1):
					for k in xrange(j+1):
						self.idx[0]=self.idxs[n,0]=i
						self.idx[1]=self.idxs[n,1]=j
						self.idx[2]=self.idxs[n,2]=k	
						m=self.pa(self.idx)				
						self.Ns[m]=n
						n+=1
		if nh==3:
			for i in xrange(natom):
				for j in xrange(i+1):
					for k in xrange(j+1):
						for p in xrange(k+1):
							self.idx[0]=self.idxs[n,0]=i
							self.idx[1]=self.idxs[n,1]=j
							self.idx[2]=self.idxs[n,2]=k	
							self.idx[3]=self.idxs[n,3]=p	
							m=self.pa(self.idx)				
							self.Ns[m]=n
							n+=1			
		self.N=n
	cdef int pa(self,int[:] idx):
		cdef int s=0
		cdef int i
		nh=self.nh
		for i in xrange(nh,-1,-1):
			s+=idx[nh-i]*self.natompow[i]
		return self.Ns[s]	
	cdef int map(self,int ii,int[:] map1):
		cdef int[:] c
		cdef int jj
		c=self.idxs[ii]
		for u in xrange(len(c)):
			self.idx1[u]=map1[c[u]]
		self.idx1=self.idx1.sort()
		jj=self.pa(self.idx1)
		return jj