#encoding : utf8
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
import numpy as np
from scipy import stats
import math,sys,os
from aces.tempAve import drawTempAve
from aces.input import postMini;
from aces.fixAveSpace import fixAveSpace
import pandas as pd
class profile:
	def __init__(self):
		self.method='muller'
		pass
	def getTempProfile(self,begin,upP,deta,S,tcfactor,zfactor):
		fas=fixAveSpace('tempProfile.txt')
		quants=fas.getIbin[12].cumsum(axis=0)
		for i in range(len(quants)):
			quants[i]/=i+1
		quants=pd.DataFrame(quants,names=['Temperature(K)','jx'])
		quants.to_csv('convergenceT.txt',sep='\t',index=None)
		self.fpro=open('tempProfile.txt')
		snapStep,nbin=fas.getNbin()
		sumTemp=np.zeros(nbin)
		sumjx=np.zeros(nbin)
		sumN=np.zeros(nbin)
		n=0
		fk=open('convergenceK.txt','w')
		pl.figure()

		pl.xlabel('x(Augstrom)')
		pl.ylabel('temperature(K)')

		for istep in range(begin,fas.nstep):
			coord,ncount,quants=fas.getIStep(istep)
			v_temp,jx=quants
			n+=1
			sumTemp+=v_temp
			sumjx+=jx
			sumN+=ncount

			
			#kappa convergence
			filter=sumN>0
			aveC=coord[filter]
			aveN=sumN[filter]/float(n)
			aveTemp=sumTemp[filter]/n
			avejx=sumjx[filter]/n
			if istep%10==0:
				pl.plot(aveC,aveTemp,label="time=%s"%(snapStep*istep))
			slope,flux_bulk=self.sslope(aveC,aveTemp,aveN,avejx,upP,deta,S)
			kappa=self.getFx(istep)/slope*tcfactor*zfactor
			fk.write("%d\t%f\n"%(istep,kappa))
		#pl.legend()
		pl.savefig('convergenceT.png',bbox_inches='tight',transparent=True) 	
		filter=sumN>0
		aveC=coord[filter]
		aveN=sumN[filter]/float(n)
		aveTemp=sumTemp[filter]/n
		avejx=sumjx[filter]/n
		nbin=len(avejx)
		fave=open('tempAve.txt','w')	
		s="id\tCoord\tCount\tTemp\tJx\n"
		for i in range(nbin):
			s+="%d\t%f\t%f\t%f\t%f\n"%(i+1,aveC[i],aveN[i],aveTemp[i],avejx[i])
		fave.writelines(s)
		fave.close()
		drawTempAve()
		return (aveC,aveN,aveTemp,avejx)
	
	def nvtSlope(self,aveC,aveTemp,aveN,avejx,upP):
		m=len(aveC);
		pt1=upP;
		pt2=m-upP-1;
		savejx=avejx[pt1:pt2+1]
		saveN=aveN[pt1:pt2+1]
		ave_jx=np.average(np.abs(savejx))
		ave_N=np.average(saveN)
		J_bulk=ave_jx*ave_N
		J_bulkc=np.average(np.abs(saveN*savejx))
		slope=np.abs(self.slope(aveC,aveTemp,pt1,pt2))
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

	def mullerSlope(self,aveC,aveTemp,aveN,avejx,upP):
		m=len(aveC);
		downP=upP;
		if m%2==0:
			cter2=m/2
			cter1=cter2-1
			pt11=upP;pt12=cter1-upP
			pt21=cter2+upP;pt21=m-1-upP
		else:
			cter=int((m-1)/2)
			pt11=upP;pt12=cter-upP
			pt21=cter+upP;pt21=m-1-upP
		slope1=self.slope(aveC,aveTemp,pt11,pt12);
		savejx=avejx[pt11:pt12+1]
		saveN=aveN[pt11:pt12+1]
		ave_jx=np.average(np.abs(savejx));
		ave_N=np.average(saveN);
		J_bulk1=ave_jx*ave_N;
		J_bulkc1=np.average(np.abs(saveN*savejx));
		slope2=-self.slope(aveC,aveTemp,pt21,pt22);
		savejx=avejx[pt21:pt22+1]
		saveN=aveN[pt21:pt22+1]
		ave_jx=np.average(np.abs(savejx));
		ave_N=np.average(saveN);
		J_bulk2=ave_jx*ave_N;
		J_bulkc2=np.average(np.abs(saveN*savejx));
		slope=(slope1+slope2)/2;
		J_bulk=(J_bulk1+J_bulk2)/2;
		J_bulkc=(J_bulkc1+J_bulkc2)/2;
		return (slope,J_bulk,J_bulkc);
		
	def getFlux(self,begin,timestep,S,conti,lz,excRate,swapEnergyRate):
		fx=[]
		st=begin;
		method=self.method
		step=[]
		hot=[]
		if(method=="nvt"):
			nvtWork=open("nvtWork.txt");
			nvtWork.next()
			line=nvtWork.next();
			co=0;
			for line in nvtWork:
				f_step,f_hot=line.strip().split()[:2]
				step.append(float(f_step))
				hot.append(float(f_hot))
				if co<=st:hotslope=0
				else:
					hotslope=np.abs(hot[co]-hot[st])/(step[co]-step[st])
				J=hotslope/timestep;
				flux_src=J/S;
				fx.append(flux_src)
				co+=1

			nvtWork.close()
		
		

		if(method=="muller"):
			file=open("swapEnergy.txt");
			file.next()
			file.next()
			co=0;
			sum=0;

			for line in file:
				f_step,heat_swap=line.strip().split()[:2]
				step.append(float(f_step))
				hot.append(float(heat_swap))
				sum+=heat_swap;
				if co<=st:hotslope=0
				else:
					ave_heat_swap=np.abs(hot[co]-hot[st])/(step[co]-step[st])
				if(conti):ave_heat_swap=sum*lz/co/excRate;
				J=ave_heat_swap/(timestep);
				flux_src=J/(2*S);
				fx.append(flux_src)
				co+=1
			
			file.close();
		
		
		if(method=="inject"):
			J=swapEnergyRate;
			flux_src=J/(2*S);
			fx.append(flux_src)
			
		self.fx=fx
		self.flux_src=flux_src
		return (flux_src,fx);
	
	def getFx(self,istep):
		fx=self.fx
		if istep<len(fx):return fx[istep]
		else: 
			return fx[-1]

	def slope(self,x,y,pt1,pt2):
		n=pt2-pt1+1;
		sxy=0.0;sx=0.0;sy=0.0;sx2=0.0;
		for i in range(pt1,pt2+1):
			sxy+=x[i]*y[i];
			sx+=x[i];
			sy+=y[i];
			sx2+=x[i]*x[i];
		return (n*sxy-sx*sy)/(n*sx2-sx*sx)
	
	
		
	
def run(method,begin,timestep,conti,excRate,swapEnergyRate,upP,deta,tcfactor,fourierTc ,computeTc ,corRate ,kb ,T,xp,yp,zp,enforceThick,thick):
	#lx,ly,S,zfactor from postMini
	p=profile()
	p.method=method
	tcfactor,deta,timestep,excRate,swapEnergyRate,corRate,kb,T,thick=map(float,[tcfactor,deta,timestep,excRate,swapEnergyRate,corRate,kb,T,thick])
	upP,begin,conti,fourierTc,computeTc,xp,yp,zp,enforceThick=map(int,[upP,begin,conti,fourierTc,computeTc,xp,yp,zp,enforceThick])

	lx,ly,lz,zfactor,S=postMini(xp,yp,zp,enforceThick,thick)[:5]
	if method=='greenkubo':
		f=open('result.txt','w')
		if(computeTc):
			os.popen("tail -2000 kappa.txt>tailKp.txt 2>err");
			file=open("tailKp.txt","r");
			s=0;n=0;
			for line in file:
				step,kp=line.split()
				kp=float(kp)
				s+=kp;
				n+=1
			
			kx=s/n;
	
		elif(fourierTc):
			v=lx*ly*lz;
			factor=corRate*timestep/(v*kb*T*T)*zfactor*tcfactor;
			path=os.path.dirname(os.path.realpath(__file__))
			kx=os.popen("%s/correlation/corr factor"%path).read();#generate a correlation file named jcor.txt
		else:
		
			gk_result=os.popen("tail -1 fileKappa 2>err").read()
			kx=gk_result.split()[1]
		
		f.write("kappa_src=%f\n"%kx);
		os.exit()
	flux_src=p.getFlux(begin,timestep,S,conti,lz,excRate,swapEnergyRate)[0]

	aveC,aveN,aveTemp,avejx=p.getTempProfile(begin,upP,deta,S,tcfactor,zfactor)
	slope,flux_bulk=p.sslope(aveC,aveTemp,aveN,avejx,upP,deta,S);
	kappa_src=flux_src/slope*tcfactor*zfactor;
	kappa_bulk=flux_bulk/slope*tcfactor*zfactor;
	f=open('result.txt','w')
	f.write('method:'+method+'\n');
	f.write("kappa_src=%f\n"%(kappa_src));
	f.close()


	fileScan=open("scan.txt","w");
	fileScan.write("method:%s\n"%method);
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
	
	fileScan.write("upP\tkappa_src\tkappa_bulk\tkappa_bulkc\tflux_src\tflux_bulk\tflux_bulkc\tslope\n");
	for i in range(0,numS):
		kappa_src=flux_src/slopes[i]*tcfactor*zfactor;
		flux_bulk=J_bulks[i]/(deta*S);
		flux_bulkc=J_bulkcs[i]/(deta*S);
		kappa_bulk=flux_bulk/slopes[i]*tcfactor*zfactor;
		kappa_bulkc=flux_bulkc/slopes[i]*tcfactor*zfactor;
		fileScan.write("%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"%(i+1,kappa_src,kappa_bulk,kappa_bulkc,flux_src,flux_bulk,flux_bulkc,slopes[i]));
	
	fileScan.close()
from os.path import *	
def proc():

	import sys
	import os,json,imp
	from aces.Units import Units
	home=dirname(realpath(__file__))
	#app home 
	projHome=dirname(realpath(sys.argv[0]))
	f=open('app.json')
	opt=f.read()
	opt=json.loads(opt)
	f.close()
	species=opt['species']
	m= imp.load_source('structure', home+'/materials/'+species+'/structure.py') 
	m=m.structure(home,opt)
	units=Units(m.units)
	m.kb=units.boltz
	m.nktv=units.nktv2p
	if(m.method=="nvt"):m.xp=0;
	lx,ly,lz,m.zfactor,m.S,xlo,xhi,ylo,yhi,zlo,zhi=postMini(m.xp,m.yp,m.zp,m.enforceThick,m.thick)
	m.dtime=m.timestep*100;
	m.tcfactor=units.tcfactor;
	m.excNum=m.aveRate/m.excRate;
	m.swapEnergyRate=m.swapEnergy/(m.excRate*m.timestep);
	run(m.method,m.begin,m.timestep,m.conti,m.excRate,m.swapEnergyRate,m.upP,m.deta,m.tcfactor,m.fourierTc ,
	m.computeTc ,m.corRate ,m.kb ,m.T,m.xp,m.yp,m.zp,m.enforceThick,m.thick)	
