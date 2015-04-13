#encoding : utf8
import matplotlib
matplotlib.use('Agg')
import numpy as np
from scipy import stats
import math,sys,os
from aces.tempAve import drawTempAve
from aces.input import postMini;
class profile:
	def __init__(self):
		self.method='muller'
		pass
	def getTempProfile(self,begin,upP,deta,S,tcfactor,zfactor):
		self.fpro=open('tempProfile.txt')
		self.nbin=nbin=self.getNbin()
		sumTemp=np.zeros([nbin,1])
		sumjx=np.zeros([nbin,1])
		sumN=np.zeros([nbin,1])
		istep=-1
		n=0
		ft=open('convergenceT.txt','w')	
		ft.write("step\ttemperature\tjx\n")
		fk=open('convergenceK.txt','w')
		while self.fpro.readline():
			istep+=1
			coord,ncount,v_temp,jx=self.getBinInfo()
			if istep<begin:continue
			n+=1
			sumTemp+=v_temp
			sumjx+=jx
			sumN+=ncount
			att=sumTemp[12]/n
			atj=sumjx[12]/n
			ft.write("%d\t%f\t%f\n"%(istep,att,atj))
			
			#kappa convergence
			aveTemp1=sumTemp/n
			avejx1=sumjx/n
			aveN1=sumN/n
			slope,flux_bulk=self.sslope(coord,aveTemp1,aveN1,avejx1,upP,deta,S)
			kappa=self.getFx(istep)/slope*tcfactor*zfactor
			fk.write("%d\t%f\n"%(istep,kappa))
			
		sumTemp/=n
		sumjx/=n
		sumN/=n
		filter=sumN>0
		aveN=sumN[filter]
		aveTemp=sumTemp[filter]
		avejx=sumjx[filter]
		fave=open('tempAve.txt','w')	
		s="id\tCoord\tCount\tTemp\tJx\n"
		for i in range(nbin):
			s+="%d\t%f\t%f\t%f\t%f\n"%(i+1,coord[i],aveN[i],aveTemp[i],avejx[i])
		fave.writelines(s)
		return (coord,aveN,aveTemp,avejx)
	
	def nvtSlope(self,aveC,aveTemp,aveN,avejx,upP):
		m=len(aveC);
		downP=upP;
		pt1=downP;
		pt2=m+1-upP;
		pt1-=1
		pt2-=1
		n=pt2-pt1+1;
		savejx=avejx[pt1:pt1+n+1]
		saveN=aveN[pt1:pt1+n+1]
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
		cter=int((m+1)/2.0);
		pt11=downP;pt12=cter-upP;
		pt22=m+1-downP;pt21=cter+upP;
		pt11-=1;pt12-=1;pt21-=1;pt22-=1;
		slope1=self.slope(aveC,aveTemp,pt11,pt12);
		n=pt12-pt11+1;
		savejx=avejx[pt11:pt11+n]
		saveN=aveN[pt11:pt11+n]
		ave_jx=np.average(np.abs(savejx));
		ave_N=np.average(saveN);
		J_bulk1=ave_jx*ave_N;
		J_bulkc1=np.average(np.abs(saveN*savejx));
		slope2=-self.slope(aveC,aveTemp,pt21,pt22);
		n=pt22-pt21+1;
		savejx=avejx[pt21:pt21+n]
		saveN=aveN[pt21:pt21+n]
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
		sxy=0;sx=0;sy=0;sx2=0;
		for i in range(pt1,pt2+1):
			sxy+=x[i]*y[i];
			sx+=x[i];
			sy+=y[i];
			sx2+=x[i]*x[i];
		
		return (n*sxy-sx*sy)/(n*sx2-sx*sx)
	
	def getNbin(self):
		fpro=self.fpro
		for i in range(3):fpro.readline()
		start=fpro.tell()
		line=fpro.readline()
		timestep,nbin=line.strip().split()
		fpro.seek(start)
		return int(nbin)	
		
	def getBinInfo(self):
		fpro=self.fpro
		nbin=self.nbin
		bin=np.zeros([nbin,1])
		coord=np.zeros([nbin,1])
		ncount=np.zeros([nbin,1])
		v_temp=np.zeros([nbin,1])
		jx=np.zeros([nbin,1])
		for i in range(nbin):
			line=fpro.readline()
			bin[i],coord[i],ncount[i],v_temp[i],jx[i]=line.strip().split()
		return (coord,ncount,v_temp,jx)	
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
	n=int(lx/2/deta)+1
	rg=range(1,n)
	slopes=np.zeros([n,1])
	J_bulks=np.zeros([n,1])
	J_bulkcs=np.zeros([n,1])
	if(method=="muller" or method=="inject"):
		for upP in range(1,n/2):
			lslopes[numS],J_bulks[numS],J_bulkcs[numS]=p.mullerSlope(aveC,aveTemp,aveN,avejx,upP);
			numS+=1
		
	
	if(method=="nvt"):
		for upP in range(1,n):
			slopes[numS],J_bulks[numS],J_bulkcs[numS]=p.nvtSlope(aveC,aveTemp,aveN,avejx,upP);
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

	
	drawTempAve()		
if __name__=='__main__':

	method,begin,timestep,conti,excRate,swapEnergyRate,upP,deta,tcfactor,fourierTc ,computeTc ,corRate ,kb ,T,xp,yp,zp,enforceThick,thick=sys.argv[1:]
	run(method,begin,timestep,conti,excRate,swapEnergyRate,upP,deta,tcfactor,fourierTc ,computeTc ,corRate ,kb ,T,xp,yp,zp,enforceThick,thick)