import numpy as np
from aces.lineManager import  lineManager
class fixAveTime:
	def __init__(self,filename):
		lm=lineManager(filename)
		self.title=lm.getLine(0).replace("# ","")
		s=lm.getLine(1).replace("# ","")
		labels=s.split(' ')
		quantsLabel=labels[1:]
		nquants=len(quantsLabel)
		self.lm=lm
		nstep=lm.nline-2
		steps=np.zeros(nstep)
		quants=np.zeros([nstep,nquants])
		for istep in range(nstep):
			line=lm.nextLine().split()
			steps[istep]=line[0]
			quants[istep]=line[1:]
		self.steps=steps
		self.quants=quants
		self.nquants=nquants
		
	def getSlopes(self,begin):
		steps=self.steps
		nstep=len(steps)
		nquants=self.nquants
		quants=self.quants
		slopes=np.zeros([nstep,nquants])
		for istep in range(begin+1,nstep):
			slopes[istep]=(quants[istep]-quants[begin])/(steps[istep]-steps[begin])
		return slopes	
	def getConvergence(self,begin):
		quants=self.quants[begin:].cumsum(axis=0)
		for i in range(len(quants)):
			quants[i]/=i+1
		return np.hstack([self.steps,quants])
	
