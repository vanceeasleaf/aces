from numpy.fft import rfft, irfft
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
from aces.tools import exit
class vdos:
	def __init__(self,timestep=0.0005):
		self.timestep=timestep
		self.run()
	def construct_index(self):
		self.fi  = file('velocity.txt')
		self.line_offset=[]
		offset=0
		for line in self.fi:
			self.line_offset.append(offset)
			offset+=len(line)
		self.fi.seek(0)
	def run(self):
		self.construct_index()
		self.readinfo()
		self.calculateDos()
	def getline(self,iline):
		offset=self.line_offset[iline]
		self.fi.seek(offset)
		return self.fi.next()
	def readinfo(self):
		self.natom=int(self.getline(3).split()[0])
		t1=int(self.getline(1).split()[0])
		self.line_interval=9+self.natom
		if len(self.line_offset)<self.line_interval:
			self.interval=1
		else:
			t2=int(self.getline(1+self.line_interval).split()[0])
			self.interval=t2-t1
		self.totalStep=len(self.line_offset)/self.line_interval
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
		v=np.zeros(self.totalStep)
		for i in range(self.totalStep):
			v[i]=float(self.getline(9+id+i*self.line_interval).split()[2+coord])
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
		fo1 = file('VACF.txt','w')
		fo1.write('correlation_time_ps vcaf_x vcaf_y vcaf_z vcaf_av\n')
		for i in range(self.totalStep/2):
			fo1.write("%f"%(self.timestep*i))
			for vcf in totalVcf[i]:
				fo1.write("\t%f"%(vcf))
			fo1.write("\n")
		fo1.close()

		totalDos = np.abs(rfft(totalVcf,axis=0))
		fo2 = file('VDOS.txt','w')
		fo2.write('Freq_THz vdos_x vdos_y vdos_z vdos_av\n')
		maxFreq=2.0/self.timestep
		for i in range(self.totalStep/2):
			fo2.write("%f"%(maxFreq*i/float(self.totalStep/2)))
			for j,dos in enumerate(totalDos[i]):
				fo2.write("\t%f"%dos)
			fo2.write("\n")
		fo2.close()

		print 'VACF and VDOS caculated OK'
		self.plot(totalVcf[:self.totalStep/2],totalDos)

	def select(self,x,n=1000):
		N=len(x)
		if N<n:
			return range(0,N)
		else:
			return range(0,N,N/n)
	def plot(self,totalVcf,totalDos):
		n,m=totalVcf.shape
		time=np.array(range(0,n))*self.timestep
		pl.figure()
		xx=self.select(time)
		time=time[xx]
		totalVcf=totalVcf[xx]
		pl.plot(time,totalVcf[:,0],label="vcf_x",linewidth=0.3)
		pl.plot(time,totalVcf[:,1],label="vcf_y",linewidth=0.3)
		pl.plot(time,totalVcf[:,2],label="vcf_z",linewidth=0.3)
		pl.xlabel('Correlation Time (ps)')
		pl.ylabel('Normalized Velocity Auto Correlation Function')
		pl.legend(loc='best').get_frame().set_alpha(0.0)
		pl.savefig('VACF.png',bbox_inches='tight',transparent=True) 
		n,m=totalDos.shape
		freq=np.linspace(0,1,n)*2.0/self.timestep
		xx=self.select(freq)
		pl.figure()
		pl.xlabel('Frequency (THz)')
		pl.ylabel('Phonon Density of States')
		pl.plot(freq,totalDos[:,0],label="dos_x",linewidth=0.3)
		pl.plot(freq,totalDos[:,1],label="dos_y",linewidth=0.3)
		pl.plot(freq,totalDos[:,2],label="dos_z",linewidth=0.3)
		pl.legend(loc='best').get_frame().set_alpha(0.0)
		pl.savefig('VDOS.png',bbox_inches='tight',transparent=True) 