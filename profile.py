#encoding : utf8
import numpy as np
class profile:
	def getTempProfile(self,begin,upP,deta,S,tcfactor,zfactor):
		self.fpro=open('tempProfile.txt')
		self.nbin=getNbin()
		sumTemp=np.zeros(nbin,1)
		sumjx=np.zeros(nbin,1)
		sumN=np.zeros(nbin,1)
		istep=-1	
		for line in self.fpro.readline():
			istep+=1
			coord,ncount,v_temp,jx=self.getBinInfo()
			if istep<begin:continue
			sumTemp+=v_temp
			sumjx+=jx
			sumN+=ncount
		sumTemp/=istep
		sumjx/=istep
		sumN/=istep
		filter=sumN>0
		aveN=sumN[filter]
		aveTemp=sumTemp[filter]
		avejx=sumjx[filter]
		fave=open('tempAve1.txt','w')	
		s="id\tCoord\tCount\tTemp\tJx\n"
		for i in range(nbin):
			s+="%d\t%f\t%f\t%f\t%f\n"%(i+1,coord[i],aveN[i],aveTemp[i],avejx[i])
		fave.writelines(s)
	def getNbin(self):
		fpro=self.fpro
		for i in range(3):fpro.readline()
		start=fpro.tell()
		line=fpro.readline()
		timestep,nbin=line.strip().split()
		fpro.seek(start)
		return nbin		
		
	def getBinInfo(self):
		nbin=self.nbin
		bin=np.zeros(nbin,1)
		coord=np.zeros(nbin,1)
		ncount=np.zeros(nbin,1)
		v_temp=np.zeros(nbin,1)
		jx=np.zeros(nbin,1)
		for i in range(nbin):
			line=fpro.readline()
			bin[i],coord[i],ncount[i],v_temp[i],jx[i]=line.strip().split()
		return (coord,ncount,v_temp,jx)	
			
if __name__=='__main__':
	import sys
	begin,upP,deta,S,tcfactor,zfactor=sys.argv[1:]
	profile().getTempProfile(int(begin),int(upP),float(deta),float(S),float(tcfactor),float(zfactor))