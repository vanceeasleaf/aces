from aces.lineManager import  lineManager
class velocity:
	def __init__(self,timestep=0.0005):
		self.timestep=timestep
		self.lm=lineManager('velocity.txt')
		self.db=h5py.File('velocity.h5')
	def info(self):
		lm=self.lm
		self.natom=int(lm.getLine(3).split()[0])
		t1=int(lm.getLine(1).split()[0])
		self.line_interval=9+self.natom
		if lm.nline<self.line_interval:
			self.interval=1
		else:
			t2=int(lm.getLine(1+self.line_interval).split()[0])
			self.interval=t2-t1
		self.totalStep=lm.nline/self.line_interval
		if self.totalStep%2==1:self.totalStep-=1
		print "Atom Number=",self.natom
		print "Total step=",self.totalStep
		print "interval=",self.interval
		self.timestep*=self.interval
		self.times=np.arange(self.totalStep)*self.timestep
		maxFreq=1/2.0/self.timestep
		self.freq=np.linspace(0,1,self.totalStep/2)*maxFreq
		return self.natom,self.totalStep,self.timestep,self.freq,self.times
	def atom(self,id):
		node='/velocity_atom/%d'%id
		if not node in self.db:
			self.db[node]=self.prepare(id)
		return self.db[node]
	def prepare(self,id):
		print 'prepare velocity_atom:%d'%id
		lm=self.lm
		v=np.zeros([self.totalStep,3])
		for i in range(self.totalStep):
			v[i,:]=map(float,lm.getLine(9+id+i*self.line_interval).split()[2:5])	
		return v