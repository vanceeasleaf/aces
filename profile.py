#encoding : utf8
import numpy as np
from scipy import stats
import math,sys,os
from aces.io.lammps.tempAve import drawTempAve
from aces.io.lammps.fixAveSpace import fixAveSpace
import pandas as pd
from aces import tools
from aces.io.lammps.fixAveTime import fixAveTime
#from aces.env import *
from aces.graph import series
class profile:
	def __init__(self):
		self.method='muller'
		
	def kappaConverge(self,istep,coord,aveN,aveQuants,*para):
		#kappa convergence
		snapStep,upP,deta,S,tcfactor,zfactor,dto=para
		filter=aveN[:]>1
		aveC=coord[filter]
		aveTemp=aveQuants[filter,0]
		avejx=aveQuants[filter,1]
		if istep%10==0:
			dto['plots'].append([aveC,aveTemp,"time=%s"%(snapStep*istep)])
		slope,flux_bulk=self.sslope(aveC,aveTemp,aveN,avejx,upP,deta,S)
		kappa=self.getFx(istep)/slope*tcfactor*zfactor
		dto['log']+="%d\t%f\n"%(istep,kappa)
		
	def getTempProfile(self,begin,upP,deta,S,tcfactor,zfactor):
		fas=fixAveSpace('tempProfile.txt')
		quants=fas.getConvergence(upP+1,begin)
		tools.to_txt(['Temperature(K)','Jx'],quants,'convergenceT.txt')
		snapStep=fas.snapStep

		dto=dict(log="step\tkappa\n",plots=[])
		coord,aveN,aveQuants=fas.iterate(begin,self.kappaConverge,snapStep,upP,deta,S,tcfactor,zfactor,dto)
		tools.write(dto['log'],'convergenceK.txt')
		series('x(Augstrom)','temperature(K)',dto['plots'],'convergenceT.png',linewidth=1,legend=False)
		filter=aveN[:]>1# at least 1 atom in the bin
		aveC=coord[filter]
		aveN=aveN[filter]
		aveTemp=aveQuants[filter,0]
		avejx=aveQuants[filter,1]
		nbin=len(avejx)
		data=np.c_[np.arange(nbin)+1,aveC,aveN,aveTemp,avejx]

		tools.to_txt(['id','Coord','Count','Temp','Jx'],data,'tempAve.txt')
		drawTempAve()
		return (aveC,aveN,aveTemp,avejx)
	
	def nvtSlope(self,aveC,aveTemp,aveN,avejx,upP):
		m=len(aveC);
		pt1=upP;
		pt2=m-upP-1;
		slope,J_bulk,J_bulkc=self.getSideSlope(aveC,aveTemp,aveN,avejx,pt1,pt2)
		return (slope,J_bulk,J_bulkc)
		
	def sslope(self,aveC,aveTemp,aveN,avejx,upP,deta,S):
		method=self.method
		if(method=="nvt"):
			slope,J_bulk=self.nvtSlope(aveC,aveTemp,aveN,avejx,upP)[:2]
			flux_bulk=J_bulk/(deta*S);
			
		if(method=="muller" or method=="inject"):
			slope,J_bulk=self.mullerSlope(aveC,aveTemp,aveN,avejx,upP)[:2]
			flux_bulk=J_bulk/(deta*S);
		
		return (slope,flux_bulk)
	def getSideRange(self,size,upP):
		m=size;
		if m%2==0:
			cter2=m/2
			cter1=cter2-1
			pt11=upP;pt12=cter1-upP
			pt21=cter2+upP;pt22=m-1-upP
		else:
			cter=int((m-1)/2)
			pt11=upP;pt12=cter-upP
			pt21=cter+upP;pt22=m-1-upP
		return (pt11,pt12,pt21,pt22)
	def getSideSlope(self,aveC,aveTemp,aveN,avejx,pt11,pt12):
		slope1=self.slopeDiff(aveC,aveTemp,pt11,pt12)
		slope1=np.abs(slope1)
		savejx=avejx[pt11:pt12+1]
		saveN=aveN[pt11:pt12+1]
		ave_jx=np.average(np.abs(savejx));
		ave_N=np.average(saveN);
		J_bulk1=ave_jx*ave_N;
		J_bulkc1=np.average(np.abs(saveN*savejx))
		return (slope1,J_bulk1,J_bulkc1)
	def mullerSlope(self,aveC,aveTemp,aveN,avejx,upP):
		pt11,pt12,pt21,pt22=self.getSideRange(len(aveC),upP)
		slope1,J_bulk1,J_bulkc1=self.getSideSlope(aveC,aveTemp,aveN,avejx,pt11,pt12)
		slope2,J_bulk2,J_bulkc2=self.getSideSlope(aveC,aveTemp,aveN,avejx,pt21,pt22)
		slope=(slope1+slope2)/2;
		J_bulk=(J_bulk1+J_bulk2)/2;
		J_bulkc=(J_bulkc1+J_bulkc2)/2;
		return (slope,J_bulk,J_bulkc);
		
	def getFlux(self,begin,timestep,S,conti,lz,excRate,swapEnergyRate):
		method=self.method
		if(method=="nvt"):
			fat=fixAveTime("nvtWork.txt")
			fx=.5*(fat.getSlopes(begin)[:,0]-fat.getSlopes(begin)[:,1])
			fx=np.abs(fx)/timestep/S
			flux_src=fx[-1]
		

		if(method=="muller"):
			fat=fixAveTime("swap.txt")
			fx=fat.getSlopes(begin)[:,0]
			fx=np.abs(fx)/timestep/S/2
			
			if(conti):
				fx=fat.getConvergence(begin)[:,0]
				fx=np.abs(fx)*lz/excRate/timestep/S/2
			flux_src=fx[-1]	
			
		if(method=="inject"):
			J=swapEnergyRate;
			flux_src=J/(2*S);
			fx=[flux_src]
			
		self.fx=fx
		return (flux_src,fx);
	
	def getFx(self,istep):
		fx=self.fx
		if istep<len(fx):return fx[istep]
		else: 
			return fx[-1]
	def slopeDiff(self,x,y,pt1,pt2):
		return (y[pt2]-y[pt1])/(x[pt2]-x[pt1])
	def slope(self,x,y,pt1,pt2):
		n=pt2-pt1+1;
		sxy=0.0;sx=0.0;sy=0.0;sx2=0.0;
		for i in range(pt1,pt2+1):
			sxy+=x[i]*y[i];
			sx+=x[i];
			sy+=y[i];
			sx2+=x[i]*x[i];
		return (n*sxy-sx*sy)/(n*sx2-sx*sx)
	
	
		
	
def run(method,begin,timestep,conti,excRate,swapEnergyRate,upP,deta,tcfactor,fourierTc ,computeTc ,corRate ,kb ,T,xp,yp,zp,enforceThick,thick,zfactor,S,box,**rest):

	p=profile()
	p.method=method
	lx,ly,lz=box[-3:]
	if method=='greenkubo':
		f=open('result.txt','w')
		if(computeTc):
			os.popen("tail -2000 kappa.txt>tailKp.txt 2>err");
			df=pd.read_csv("tailKp.txt",sep=' ',header=None)
			kx=np.average(np.array(df)[:,1],axis=0)
	
		elif(fourierTc):
			v=lx*ly*lz;
			factor=corRate*timestep/(v*kb*T*T)*zfactor*tcfactor;
			SRCHOME=dirname(__FILE__)
			kx=os.popen("%s/correlation/corr factor"%SRCHOME).read();#generate a correlation file named jcor.txt
		else:
		
			gk_result=os.popen("tail -1 fileKappa 2>err").read()
			kx=gk_result.split()[1]
		
		f.write("kappa_src=%f\n"%kx);
		return
	flux_src=p.getFlux(begin,timestep,S,conti,lz,excRate,swapEnergyRate)[0]

	aveC,aveN,aveTemp,avejx=p.getTempProfile(begin,upP,deta,S,tcfactor,zfactor)
	slope,flux_bulk=p.sslope(aveC,aveTemp,aveN,avejx,upP,deta,S);
	kappa_src=flux_src/slope*tcfactor*zfactor;
	kappa_bulk=flux_bulk/slope*tcfactor*zfactor;
	f=open('result.txt','w')
	f.write('method:'+method+'\n');
	f.write("kappa_src=%f\n"%(kappa_src));
	f.close()



	numS=0;
	n=len(aveC)-3
	slopes=np.zeros(n)
	J_bulks=np.zeros(n)
	J_bulkcs=np.zeros(n)
	if(method=="muller" or method=="inject"):
		for upP in range(1,n/4):
			s=p.mullerSlope(aveC,aveTemp,aveN,avejx,upP);
			slopes[numS],J_bulks[numS],J_bulkcs[numS]=s
			numS+=1

		
	
	if(method=="nvt"):
		for upP in range(1,n/2):
			s=p.nvtSlope(aveC,aveTemp,aveN,avejx,upP);
			slopes[numS],J_bulks[numS],J_bulkcs[numS]=s
			numS+=1
	

	kappa_src=flux_src/slopes*tcfactor*zfactor;
	flux_bulk=J_bulks/(deta*S);
	flux_bulkc=J_bulkcs/(deta*S);
	kappa_bulk=flux_bulk/slopes*tcfactor*zfactor;
	kappa_bulkc=flux_bulkc/slopes*tcfactor*zfactor;
	data=np.c_[np.arange(n)+1,kappa_src,kappa_bulk,kappa_bulkc,np.ones(n)*flux_src,flux_bulk,flux_bulkc,slopes]
	tools.to_txt('upP kappa_src kappa_bulk kappa_bulkc flux_src flux_bulk flux_bulkc slope'.split(' '),data[:numS],"scan.txt")
	

