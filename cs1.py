from aces.cs import *
class runner:
	def __init__(self,NAH=2,split=True,mu=0.1,lam=0.9):
		self.NAH=NAH
		self.split=split
		self.mu=mu
		self.lam=lam
	def getTrainSets(self,u):
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

		F=np.load('allforce.npy')
		u=np.load('allpos.npy')-io.read('POSCAR').positions[np.newaxis,:,:]
		v=norm(u,axis=2)
		u0=v.flatten().max()

		F,A=self.getMatrix(F,u/u0)
		print "start compressive sensing "
		B=cssklearn().run(F,A)
		assert not np.allclose(B,0)
		print "rebuilding IFCs "
		phi=self.rebuild(B)
		print "writing IFCs "
		fc2=np.einsum(phi[1],[1,0,3,2])/u0

		writefc2(fc2,'csfc2')
	def rebuild(self,B):
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