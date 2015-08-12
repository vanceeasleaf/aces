from numpy.fft import rfft, irfft
import numpy as np
from aces.graph import series,plot
from aces.tools import exit,write,to_txt
from aces.lineManager import  lineManager
from scipy.optimize import leastsq
from aces.dos import plot_smooth,plot_dos,plot_vacf,plot_atomdos
from math import pi
import h5py
class vdos:
	def __init__(self,timestep=0.0005):
		self.timestep=timestep
		self.phase=None
		print "scanning the velocity file"
		self.lm=lineManager('velocity.txt')
		#self.run()
		self.db=h5py.File('velocity.hdf5')
		self.readinfo()
	def run(self):

		
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
		self.times=np.arange(self.totalStep)*self.timestep
		maxFreq=1/2.0/self.timestep
		self.freq=np.linspace(0,1,self.totalStep/2)*maxFreq
	def correlate(self,a,b):
		length = len(a)
		if a==b:
			b = rfft(b,axis=0)
			a = b
		else:
			b = rfft(b,axis=0)
			a=rfft(a,axis=0)
		a = a.conjugate()     #  a(t0)b(t0+t)
		c = irfft(a*b,axis=0)/length
		return c
	def correlation_atom(self,id):
		node='/correlate_atom/%d'%id
		if not node in self.db:
			print 'prepare correlate_atom:%d'%id

			v=self.velocity_atom(id)
			self.db[node]=self.correlate(v,v)
		return self.db[node]
	def dos_atom(self,id):
		if not '/freq' in self.db:	
			self.db['/freq']= self.freq
		node='/dos_atom/%d'%id
		if not node in self.db:
			print 'prepare dos_atom:%d'%id
			vcf=self.correlation_atom(id)
			vcf=np.average(vcf[:,:3],axis=1)
			vcf0=vcf[0].copy()
			vcf/=vcf0			
			dos=np.abs(rfft(vcf,axis=0))			
			self.db[node]= dos[:self.totalStep/2]
		return self.db[node]
	def cal_atomdos(self):
		for i in range(self.natom):
			self.dos_atom(i)
		plot_atomdos()
	def calculateDos(self):
		totalVcf=np.zeros([self.totalStep,4])
		for i in range(self.natom):
			print "atom",i
			vcf=self.correlation_atom(i)
			totalVcf[:,:3]+=vcf
			totalVcf[:,3]=np.average(totalVcf[:,:3],axis=1)
		#normalization
		vcf0=totalVcf[0].copy()
		for i,vcf in enumerate(totalVcf):
			totalVcf[i]/=vcf0

		totalStep=self.totalStep
		data=np.c_[self.times,totalVcf]
		to_txt(['correlation_time(ps)','vcaf_x','vcaf_y','vcaf_z','vcaf_av'],data[:totalStep/2],'VACF.txt')
		totalDos = np.abs(rfft(totalVcf,axis=0))[:totalStep/2]
		

		data=np.c_[self.freq,totalDos]
		to_txt(['Freq_THz','vdos_x','vdos_y','vdos_z','vdos_av'],data,'VDOS.txt')
		
		print 'VACF and VDOS caculated OK'
		plot_dos()
		plot_vacf()
		plot_smooth()
	def velocity_atom(self,id):
		node='/velocity_atom/%d'%id
		if not node in self.db:
			print 'prepare velocity_atom:%d'%id
			lm=self.lm
			v=np.zeros([self.totalStep,3])
			for i in range(self.totalStep):
				v[i,:]=map(float,lm.getLine(9+id+i*self.line_interval).split()[2:5])	
			self.db[node]=v
		return self.db[node]
	def fourier_atom(self,id):
		node='/fourier_atom/%d'%id
		if not node in self.db:
			print 'prepare fourier_atom:%d'%id
			v=self.velocity_atom(id)
			r=rfft(v,axis=0)
			self.db[node]=r
		return self.db[node]
	def getPhase(self,k,natom_unitcell):
		p=np.exp(self.pfactor.dot(k))
		v=p[range(0,self.natom,natom_unitcell)]
		self.phase=np.repeat(v,natom_unitcell)
		return self.phase

	def calculateLife(self,eigen):
		totalStep=self.totalStep
		k,freq,vec=eigen
		vec=vec.conjugate()
		q=np.zeros(totalStep/2+1,dtype=np.complex)
		nu=len(vec)
		phase=self.getPhase(k,nu)
		for i in range(self.natom):
			
			fv=self.fourier_atom(i)
			iu=i%nu
			q+=np.array(fv).dot(vec[iu])*phase[i]
		q=(q*q.conjugate()).real
		
		x=np.linspace(0,1,totalStep/2+1)*1/2.0/self.timestep
		low=max(q.argmax()-20,0)
		hi=min(q.argmax()+20,len(q))
		filter=range(low,hi)
		x=x[filter]
		q=q[filter]
		p0 = np.array([x[q.argmax()], 0.1, q[q.argmax()]]) #Initial guess
		p=self.fitLife(x,q,p0)
		#to_txt(['Freq','dos'],np.c_[x,q],'single%s%s.txt'%(str(k),freq))
		"""
		series(xlabel='Frequency (THz)',
			ylabel='Single Phonon Power Spectrum',
			datas=[(x,q,"origin"),
			(x,self.lorentz(p,x),"fitting")]
			,linewidth=1
			,filename='single%s%s.png'%(str(k),freq))
		"""
		v=map(str,list(k)+[freq]+list(p))
		return '\t'.join(v)
	def life_yaml(self,filename="mesh.yaml",correlation_supercell=[10,10,1]):
		from aces.lammpsdata import lammpsdata
		atoms=lammpsdata().set_src('correlation_structure')
		self.pfactor=1j*atoms.positions.dot(atoms.get_reciprocal_cell()*np.c_[correlation_supercell])*2*pi
		
		from aces.phononyaml import phononyaml
		
		pya=phononyaml(filename)
		f=open('life.txt','w')
		f.write('kx\tky\tkz\tfreq\tw0\ttao\ta0\n')
		for iqp in range(pya.nqpoint):
			for ibr in range(pya.nbranch):
				print 'branch:%s %s'%(iqp,ibr)
				
				k=pya.qposition(iqp)
				freq=pya.frequency(iqp,ibr)
				vec=pya.atoms(iqp,ibr)
				v=self.calculateLife((k,freq,vec))
				f.write('%s\n'%v)


	def fitLife(self,x,z,p0):
		

		solp, ier = leastsq(self.errorfunc, 
                    p0, 
                    args=(x,z),
                    Dfun=None,
                    full_output=False,
                    ftol=1e-9,
                    xtol=1e-9,
                    maxfev=100000,
                    epsfcn=1e-10,
                    factor=0.1)
		return solp
	def lorentz(self,p,x):
		return p[2] / ((x-p[0])**2 + p[1]**2/4)

	def errorfunc(self,p,x,z):
		return self.lorentz(p,x)-z	
		
