from numpy.fft import rfft, irfft
import numpy as np
from aces.graph import series,plot
from aces.tools import exit,write,to_txt
from aces.velocityh5 import velocity
from scipy.optimize import leastsq
from aces.dos import plot_smooth,plot_dos,plot_vacf,plot_atomdos
from math import pi
from aces.tools import *
from scipy import signal
import h5py

from aces.qpointsyaml import phononyaml
import time
class vdos:
	def __init__(self,timestep=0.0005):
		self.phase=None
		self.velocity=velocity(timestep)
		#self.run()
		self.db=h5py.File('dos.h5')
		self.dbi=0
		self.readinfo()
		self.totalsed=False
		self.partsed=False
	def run(self):

		
		self.calculateDos()

	def readinfo(self):
		self.natom,self.totalStep,self.timestep,self.freq,self.times=self.velocity.info()
	def acf(self,a):
		length = len(a)
		a = rfft(a,axis=0)
		c = irfft(a*a.conjugate(),axis=0)/length
		return c
	def correlation_atom(self,id):
		node='/correlate_atom/%d'%id
		if not node in self.db:
			print 'prepare correlate_atom:%d'%id

			v=self.velocity_atom(id)
			self.db[node]=self.acf(v)
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
		return self.velocity.atom(id)
	def fourier_atom(self,id):
		node='/fourier_atom/%d'%id
		self.dbi+=1
		if self.dbi%100000==0:
			self.db.close()
			self.db=h5py.File('dos.h5')
			self.dbi=0
			
		if not node in self.db:
			print 'prepare fourier_atom:%d'%id
			v=self.velocity_atom(id)
			r=rfft(v,axis=0)
			self.db[node]=r
		return self.db[node][:]
	def getPhase(self,k,natom_unitcell):
		p=np.exp(self.pfactor.dot(k))
		v=p[range(0,self.natom,natom_unitcell)]
		self.phase=np.repeat(v,natom_unitcell)
		return self.phase
	def calculateSED(self,k,natom_unitcell):
		totalStep=self.totalStep
		phi=np.zeros(totalStep/2+1)
		phase=self.getPhase(k,natom_unitcell)
		
		for j in range(natom_unitcell):
			q=np.zeros([totalStep/2+1,3],dtype=np.complex)
			for i in range(j,self.natom,natom_unitcell):
				fv=self.fourier_atom(i)
				q+=np.array(fv)*phase[i]
			phi+=(q*q.conjugate()).real.sum(axis=1)
		if self.partsed:
			x=np.linspace(0,1,totalStep/2+1)*1/2.0/self.timestep

			series(xlabel='Frequency (THz)',
				ylabel='Single Phonon Power Spectrum',
				datas=[(x,phi,"origin")]
				,linewidth=1
				,filename='SEDBAND/single%s.png'%(str(k)))
		return phi

	def calculateLife(self,eigen):
		totalStep=self.totalStep
		k,freq,vec=eigen
		vec=vec.conjugate()
		q=np.zeros(totalStep/2+1,dtype=np.complex)
		nu=len(vec)
		natom_unitcell=nu
		phase=self.getPhase(k,nu)
		for j in range(natom_unitcell):
			for i in range(j,self.natom,natom_unitcell):
				fv=self.fourier_atom(i)
				q+=np.array(fv).dot(vec[j])*phase[i]
		q=(q*q.conjugate()).real
		
		x=np.linspace(0,1,totalStep/2+1)*1/2.0/self.timestep
		q1=q
		x1=x
		low=max(q.argmax()-100,0)
		hi=min(q.argmax()+100,len(q))
		filter=range(low,hi)
		x=x[filter]
		q=q[filter]
		p0 = np.array([x[q.argmax()], 0.01, q[q.argmax()]]) #Initial guess
		p=self.fitLife(x,q,p0)
		#to_txt(['Freq','dos'],np.c_[x,q],'SED/single%s%s.txt'%(str(k),freq))

		if self.partsed:
			xv=np.linspace(x.min(),x.max(),100)
			series(xlabel='Frequency (THz)',
				ylabel='Single Phonon Power Spectrum',
				datas=[(x,q,"origin"),
				(xv,self.lorentz(p,xv),"fitting")]
				,linewidth=1
				,filename='SED/part%s%s.png'%(str(k),freq))
		if self.totalsed:
			series(xlabel='Frequency (THz)',
				ylabel='Single Phonon Power Spectrum',
				datas=[(x1,q1,"origin"),
				(xv,self.lorentz(p,xv),"fitting")]
				,linewidth=1
				,filename='SED/all%s%s.png'%(str(k),freq))
		
		v=map(str,list(k)+[freq,p[0],p[1],1.0/p[1] ])
		return '\t'.join(v)
	def getpfactor(self,correlation_supercell=[10,10,1]):
		from aces.lammpsdata import lammpsdata
		atoms=lammpsdata().set_src('correlation_structure')
		b=atoms.get_reciprocal_cell()*np.c_[correlation_supercell]
		return 1j*atoms.positions.dot(b.T)*2*pi
	def specialk(self,pya,k0=[-0.0327869,0,0]):
		self.partsed=True;
		self.totalsed=True;
		
		iqp=0
		for iqp0 in range(pya.nqpoint):
			k=pya.qposition(iqp0)
			iqp=iqp0
			if np.allclose(k0,k):
				break
		else:
			print "no k found"
			return 0
		q=self.lifeSED(k,pya.natom,iqp,pya)

		for ibr in range(pya.nbranch):
			freq=pya.frequency(iqp,ibr)
			vec=pya.atoms(iqp,ibr)
			v=self.calculateLife((k,freq,vec))
			print v
		return 0
	def life_yaml(self,filename="qpoints/qpoints.yaml",correlation_supercell=[10,10,1]):
		self.pfactor=self.getpfactor(correlation_supercell)

		if not exists('SED'):mkdir('SED')
		pya=phononyaml(filename)
		if False:
			self.specialk(pya,[-0.0327869,0,0])
			self.specialk(pya,[0.0327869,0,0])
			return 
			self.specialk(pya,[0.1333333,0,0])
			self.specialk(pya,[-0.1333333,0,0])
			return
		self.specialk(pya,[0.125,0.75,0])

		return
		f=open('life.txt','w')
		c=h5py.File('life.h5')
		f.write('kx\tky\tkz\tfreq\tw0\tscatter\ttao\n')
		for iqp in range(pya.nqpoint):
			for ibr in range(pya.nbranch):
				node='/%s/%s'%(iqp,ibr)
				print node
				if not node in c:
					k=pya.qposition(iqp)
					freq=pya.frequency(iqp,ibr)
					vec=pya.atoms(iqp,ibr)
					v=self.calculateLife((k,freq,vec))
					c[node]=v
				f.write('%s\n'%c[node][()])
	def lifesed(self,filename="qpoints/qpoints.yaml",correlation_supercell=[10,10,1]):

		self.pfactor=self.getpfactor(correlation_supercell)
		pya=phononyaml(filename)
		#self.specialk(pya,[0.5,0.0,0])

		#return
		f=open('lifesed.txt','w')
		c=h5py.File('lifesed.h5')
		f.write('kx\tky\tkz\tw\tw0\tscatter\ttao\n')
		for iqp in range(pya.nqpoint):
				node='/%s'%(iqp)
				print node
				if not node in c:
					k=pya.qposition(iqp)
					v=self.lifeSED(k,pya.natom,iqp,pya)
					c[node]=v
				f.write('%s\n'%c[node][()])
	def fitpart(self,x,q):
		p0 = np.array([x[q.argmax()], 0.01, q[q.argmax()]]) #Initial guess
		p=self.fitLife(x,q,p0)
		return p
	def lifeSED(self,k,na,iqp,pya):

		q=self.calculateSED(k,na)
		x=self.freq
		nbranch=3*na
		life=[]
		for i in range(nbranch):				
			low=max(q.argmax()-100,0)
			hi=min(q.argmax()+100,len(q)-1)
			filter=range(low,hi)
			p=self.fitpart(x[filter],q[filter])
			life.append(p)
			q[filter]=0.0
		life=np.array(life)
		w=[pya.frequency(iqp,ibr) for ibr in range(nbranch)]
		o=life[:,0].argsort()
		life=life[o]
		v=['\t'.join(map(str,list(k)+[w[i],p[0],p[1],1.0/p[1]]) ) for i,p in enumerate(life)]
		return '\n'.join(v)




	def sed_band(self,correlation_supercell=[10,10,1]):
		self.pfactor=self.getpfactor(correlation_supercell)
		data =parseyaml('band.yaml')
		natom_unitcell=int(data['natom'])
		nqpoint=int(data['nqpoint'])
		if not exists('SEDBAND'):mkdir('SEDBAND')
		sed=[]
		db=h5py.File('sed.h5')
		def u(phonon):
			qp=phonon['q-position']
			qps=str(qp)
			print qps
			if not qps in db:
				db[qps]=self.calculateSED(map(float,qp),natom_unitcell)
			s=db[qps][:]
			sed.append(s)
		map(u,data['phonon'])
		totalStep=self.totalStep
		x=np.linspace(0,1,totalStep/2+1)*1/2.0/self.timestep
		np.save('wtick.npy',x)
		np.save('sed.npy',sed)
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
		return p[2] / ((x-p[0])**2 + p[1]**(2)/4.0)

	def errorfunc(self,p,x,z):
		return self.lorentz(p,x)-z
		
