import numpy as np
from ..lineManager import  lineManager
class fixAveSpace:
	def __init__(self,filename):
		lm=lineManager(filename,cache=False)
		self.title=lm.getLine(0).replace("# ","")
		s=lm.getLine(2).replace("# ","")
		labels=s.split(' ')
		self.quants=labels[3:]
		self.nquants=len(self.quants)
		self.lm=lm
		self.snapStep,self.nbin=self.getNbin()
		self.nstep=(lm.nline-3)/(self.nbin+1)
		self.istep=0
		
	def getNbin(self):
		lm=self.lm
		line=lm.getLine(3)
		snapStep,nbin=line.split()
		return (int(snapStep),int(nbin))
		
	def nextStep(self):
		return self.getIStep(self.istep)
		
	def getIStep(self,istep):
		nstep=self.nstep
		lm=self.lm
		self.istep=istep+1
		nbin=self.nbin
		assert istep<nstep and istep>=0
		lm.getLine(3+istep*(1+nbin))

		coord=np.zeros(nbin)
		ncount=np.zeros(nbin)
		n=self.nquants
		quants=np.zeros([nbin,n])
		for i in range(nbin):
			line=lm.nextLine().split()			
			coord[i],ncount[i]=line[1:3]
			quants[i]=line[3:]
		return (coord,ncount,quants)
		
	def getIbin(self,ibin):
		lm=self.lm
		nbin=self.nbin
		nstep=self.nstep
		n=self.nquants		
		assert ibin<nbin and ibin>=0

		quants=np.zeros([nstep,n])
		for istep in range(nstep):
			line=lm.getLine(3+istep*(1+nbin)+(1+ibin)).split()
			quants[istep]=line[3:]
		return quants
		
	def getConvergence(self,ibin,begin):
		quants=self.getIbin(ibin)[begin:].cumsum(axis=0)
		for i in range(len(quants)):
			quants[i]/=i+1
		return quants
	
	def iterate(self,begin,callback,*para):
		nbin=self.nbin
		nquants=self.nquants	
		nstep=self.nstep
		sumQuants=np.zeros([nbin,nquants])
		sumN=np.zeros(nbin)
		n=0
		for istep in range(begin,nstep):
			coord,ncount,quants=self.getIStep(istep)
			n+=1
			sumQuants+=quants
			sumN+=ncount

			callback(istep,coord,sumN/n,sumQuants/n,*para)
		return (coord,sumN/n,sumQuants/n)