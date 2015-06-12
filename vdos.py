from numpy.fft import rfft, irfft
import numpy as np
from aces.graph import series
from aces.tools import exit,write,to_txt
from aces.lineManager import  lineManager
class vdos:
	def __init__(self,timestep=0.0005):
		self.timestep=timestep
		self.lm=lineManager('velocity.txt')
		self.run()

	def run(self):

		self.readinfo()
		self.calculateDos()

	def readinfo(self):
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
	def correlate(self,a,b):
		length = len(a)
		a = rfft(a).conjugate()     #  a(t0)b(t0+t)
		b = rfft(b)                 # .conjugate() for b(t0)a(t0+t)
		c = irfft(a*b)/length
		return c
	def correlation_atom(self,id,coord):
		lm=self.lm
		v=np.zeros(self.totalStep)
		for i in range(self.totalStep):
			v[i]=float(lm.getline(9+id+i*self.line_interval).split()[2+coord])
		vcf=self.correlate(v,v)
		return np.transpose(vcf)
	def calculateDos(self):
		totalVcf=np.zeros([self.totalStep,4])
		for i in range(self.natom):
			print "atom",i
			for coord in range(3):
				vcf=self.correlation_atom(i,coord)
				
				totalVcf[:,coord]+=vcf
			totalVcf[:,3]=np.average(totalVcf[:,0:4],axis=1)
		#normalization
		vcf0=totalVcf[0].copy()
		for i,vcf in enumerate(totalVcf):
			totalVcf[i]/=vcf0

		totalStep=self.totalStep
		data=np.hstack([np.arange(totalStep)*self.timestep,totalVcf])
		to_txt(['correlation_time_ps','vcaf_x','vcaf_y','vcaf_z','vcaf_av'],data[:totalStep/2],'VACF.txt')
		totalDos = np.abs(rfft(totalVcf,axis=0))
		

		maxFreq=2.0/self.timestep
		data=np.hstack([np.arange(totalStep)*maxFreq/float(self.totalStep/2),totalDos])
		to_txt(['Freq_THz','tvdos_x','tvdos_y','tvdos_z','tvdos_av'],data[:totalStep/2],'VDOS.txt')
		
		print 'VACF and VDOS caculated OK'
		self.plot(totalVcf[:totalStep/2],totalDos)

	def select(self,x,n=1000):
		N=len(x)
		if N<n:
			return range(0,N)
		else:
			return range(0,N,N/n)
	def plot(self,totalVcf,totalDos):

		n,m=totalVcf.shape
		time=np.array(range(0,n))*self.timestep
		xx=self.select(time)
		time=time[xx]
		totalVcf=totalVcf[xx]
		series(xlabel='Correlation Time (ps)',
			ylabel='Normalized Velocity Auto Correlation Function',
			datas=[(time,totalVcf[:,0],"vcf_x"),
			(time,totalVcf[:,1],"vcf_y"),
			(time,totalVcf[:,2],"vcf_z")
			,linewidth=0.3
			,filename='VACF.png')
		n,m=totalDos.shape
		freq=np.linspace(0,1,n)*2.0/self.timestep
		series(xlabel='Frequency (THz)',
			ylabel='Phonon Density of States',
			datas=[(freq,totalDos[:,0],"dos_x"),
			(freq,totalDos[:,1],"dos_y"),
			(freq,totalDos[:,2],"dos_z")
			,linewidth=0.3
			,filename='VDOS.png')