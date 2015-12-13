import h5py
import numpy as np
from aces.tools import exists,write,read
class velocity:
	def __init__(self,timestep=0.0005):
		if exists('velocity.h5md'):
			self.db=h5py.File('velocity.h5md')
	def prepare(self):
		db=self.db['/particles/main/velocity/']
		self.natom=db['value'].shape[1]
		self.totalStep=db['step'].shape[0]
		if self.totalStep%2==1:self.totalStep-=1
		print "Atom Number=",self.natom
		print "Total step=",self.totalStep
		self.timestep=db['time'][1]-db['time'][0]
		write('%d %d %f'%(self.natom,self.totalStep,self.timestep),'velocity.info')
	def info(self):
		if not exists('velocity.info'):
			self.prepare()
		s=read('velocity.info').split()
		self.natom=int(s[0])
		self.totalStep=int(s[1])
		self.timestep=float(s[2])
		self.times=np.arange(self.totalStep)*self.timestep
		maxFreq=1/2.0/self.timestep
		self.freq=np.linspace(0,1,self.totalStep/2)*maxFreq
		return self.natom,self.totalStep,self.timestep,self.freq,self.times
	def atom(self,id):
		db=self.db['/particles/main/velocity/value']
		return db[:self.totalStep,id]