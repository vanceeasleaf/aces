from numpy.fft import rfft, irfft,fft
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
		
		self.readinfo()
		self.totalsed=False
		self.partsed=False
		self.test=False
	def run(self):

		
		self.calculateDos()

	def readinfo(self):
		self.natom,self.totalStep,self.timestep,self.freq,self.times=self.velocity.info()
		self.db=h5py.File('dos.h5')
		
		
		
	def inith5(self,dataset,N,isComplex=False):
		node='/%s'%dataset
		label='/label_%s'%dataset
		if not node in self.db:
			print 'getting %s in h5'%node
			type=[np.float,np.complex][isComplex]
			self.db.create_dataset(node, (self.natom,N,3),dtype=type)
			self.db[label]=np.zeros(self.natom,dtype=np.int)# weather the atom is calculated
		return node,label
	def acf(self,a):
		length = len(a)
		a = rfft(a,axis=0)
		c = irfft(a*a.conjugate(),axis=0)/length
		return c
	def correlate_atom(self,id):
		node,label=self.inith5('correlate_atom',self.totalStep)
		if not self.db[label][id]:
			print 'prepare %s:%d'%(node,id)

			v=self.velocity_atom(id)
			self.db[node][id]=self.acf(v)
			self.db[label][id]=1
		return self.db[node][id]
	def dos_atom(self,id):
		node='/freq'
		if not node in self.db:	
			print 'getting %s in h5'%node
			self.db[node]= self.freq

		node,label=self.inith5('dos_atom',self.totalStep/2+1)
		if not self.db[label][id]:
			print 'prepare %s:%d'%(node,id)
			vcf=self.correlate_atom(id)
			vcf=np.average(vcf[:,:3],axis=1)
			vcf0=vcf[0].copy()
			vcf/=vcf0			
			dos=np.abs(rfft(vcf,axis=0))			
			self.db[node][id]= dos
			self.db[label][id]=1
		return self.db[node][id]
	def cal_atomdos(self):
		for i in range(self.natom):
			self.dos_atom(i)
		plot_atomdos()
	def calculateDos(self):
		totalVcf=np.zeros([self.totalStep,4])
		for i in range(self.natom):
			print "atom",i
			vcf=self.correlate_atom(i)
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
	def velocity_atom_harmonic(self,id,iseed):
		if self.testp is None:			
			k0=self.testk
			times=np.arange(self.totalStep)*self.timestep
			pya=self.pya
			iqp=0
			for iqp0 in range(pya.nqpoint):
				k=pya.qposition(iqp0)
				iqp=iqp0
				if np.allclose(k0,k):
					break
			else:
				print "no k found"
				return 0
			
			sphases=[]
			vecs=[]
			for ibr in range(pya.nbranch):
				freq=pya.frequency(iqp,ibr)
				vec=pya.atoms(iqp,ibr)
				natom_unitcell=len(vec)
				phase=self.getPhase(k,natom_unitcell).conjugate()
				vec=np.einsum('ij,i->ij',np.tile(vec,[self.natom/natom_unitcell,1]),phase)
				sphases.append(np.exp(1j*2.0*pi*freq*times))
				vecs.append(vec)
			self.sphases=sphases
			self.vecs=vecs
			p=np.random.rand(100,pya.nbranch)
			p=np.exp(2j*pi*p)
			self.rp=p
			self.testp=True
			self.nbranch=pya.nbranch
		p=self.rp		

		sphases=self.sphases
		vecs=self.vecs
		v=np.zeros([self.totalStep,3],dtype=np.complex)
		for ibr in range(self.nbranch):
			for i in range(3):
				vec=vecs[ibr]
				v[:,i]+=vec[id,i]*sphases[ibr]*p[iseed,ibr]
		return v.real

	def fourier_atom(self,id):
		node,label=self.inith5('fourier_atom',self.totalStep/2+1,True)
			
		if not self.db[label][id]:
			print 'prepare %s:%d'%(node,id)
			v=self.velocity_atom(id)
			r=rfft(v,axis=0)
			self.db[node][id]=r
			self.db[label][id]=1
		return self.db[node][id]
	def getPhase(self,k,natom_unitcell):
		p=np.exp(self.pfactor.dot(k))
		v=p[range(0,self.natom,natom_unitcell)]
		self.phase=np.repeat(v,natom_unitcell)
		from ase import io 
		atoms=io.read('POSCAR')
		symbols=atoms.get_chemical_symbols()
		masses=self.getMassFromLabel(symbols)
		Nc=self.natom/natom_unitcell
		masses=np.array(list(masses)*Nc)
		self.phase*=np.sqrt(masses/float(Nc))
		return self.phase
	def getMassFromLabel(self,labels):
		from ase.data import atomic_masses,atomic_numbers
		nums=[atomic_numbers[label] for label in labels]
		masses=[atomic_masses[num] for num in nums]
		return masses
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
		if self.totalsed:
			x=np.linspace(0,1,totalStep/2+1)*1/2.0/self.timestep

			series(xlabel='Frequency (THz)',
				ylabel='Phonon Energy Spectrum',
				datas=[(x,phi,"origin")]
				,linewidth=1
				,filename='NMA/sed%s.png'%(str(k)))
		return phi

	def calculateLife(self,eigen,iqp,ibr):
		totalStep=self.totalStep
		k,freq,vec=eigen
		vec=vec.real
		natom_unitcell=len(vec)
		phase=self.getPhase(k,natom_unitcell)
		vec=np.einsum('ij,i->ij',np.tile(vec,[self.natom/natom_unitcell,1]),phase)
		if not self.test:
			q=np.zeros(totalStep/2+1,dtype=np.complex)
			for i in range(self.natom):
				fv=self.fourier_atom(i)
				q+=np.array(fv).dot(vec[i])
			q=(q*q.conjugate()).real
		else:
			self.testp=None
			qc=np.zeros(totalStep/2+1)
			nseed=20
			
			for u in range(nseed):
				q=np.zeros(totalStep,dtype=np.complex)
				for i in range(self.natom):
					fv=self.velocity_atom_harmonic(i,u)
					q+=np.array(fv).dot(vec[i])
				q=fft(q)
				q=(q*q.conjugate()).real
				qc+=q[:totalStep/2+1]
			q=qc/nseed
		"""
		for j in range(natom_unitcell):
			for i in range(j,self.natom,natom_unitcell):
				fv=self.fourier_atom(i)
				q+=np.array(fv).dot(vec[j])*phase[i]
		"""
		
		x=np.linspace(0,1,totalStep/2+1)*1/2.0/self.timestep
		q1=q
		x1=x
		df=1.0/(2.0*self.timestep)/(self.totalStep/2)
		span=int(1.0/df)
		ori=int(freq/df)
		low=max(ori-1*span,0)
		hi=min(ori+1*span,len(q)-1)
		filter=range(low,hi)
		x=x[filter]
		q=q[filter]
		q=self.lowess(x,q)
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
				,filename='NMA/nmapart_%s_%s_%s_%s.png'%(iqp,ibr,str(k),freq))
		if self.totalsed:
			xv=np.linspace(x.min(),x.max(),100)
			series(xlabel='Frequency (THz)',
				ylabel='Single Phonon Power Spectrum',
				datas=[(x1,q1,"origin"),
				(xv,self.lorentz(p,xv),"fitting")]
				,linewidth=1
				,filename='NMA/nma_%s_%s_%s_%s.png'%(iqp,ibr,str(k),freq))
		
		v=map(str,list(k)+[freq,p[0],p[1],1.0/p[1] ])
		
		c='\t'.join(v)
		print "[nma]",c
		return c
	def getpfactor(self,correlation_supercell=[10,10,1]):
		from aces.lammpsdata import lammpsdata
		atoms=lammpsdata().set_src('correlation_structure')
		b=atoms.get_reciprocal_cell()*np.c_[correlation_supercell]
		return 1j*atoms.positions.dot(b.T)*2*pi
	def specialk(self,pya,k0=[0,0,0]):
		self.partsed=True
		self.totalsed=True
		iqp=0
		for iqp0 in range(pya.nqpoint):
			k=pya.qposition(iqp0)
			iqp=iqp0
			if np.allclose(k0,k):
				break
		else:
			print "no k found"
			return 0
		#q=self.lifeSED(k,pya.natom,iqp,pya)
		print "check orthorgnal"
		"""
		natom_unitcell=pya.natom
		phase=self.getPhase(k,natom_unitcell).conjugate()
		p=np.zeros([pya.nbranch]*2)
		for ibr in range(pya.nbranch):
			for ibr1 in range(pya.nbranch):
				vec=pya.atoms(iqp,ibr)
				vec1=pya.atoms(iqp,ibr1)
				vec=np.einsum('ij,i->ij',np.tile(vec,[self.natom/natom_unitcell,1]),phase)
				vec1=np.einsum('ij,i->ij',np.tile(vec1,[self.natom/natom_unitcell,1]),phase.conjugate())
				p[ibr,ibr1]=np.einsum('ij,ij',vec,vec1.conjugate())
		print p
		"""
		for ibr in range(pya.nbranch):
			freq=pya.frequency(iqp,ibr)
			vec=pya.atoms(iqp,ibr)
			v=self.calculateLife((k,freq,vec),iqp,ibr)
			
		return 0
	def lifenmaq(self,filename="qpoints/qpoints.yaml",k=[0,0,0],correlation_supercell=[10,10,1],test=False):
		self.test=test
		self.pfactor=self.getpfactor(correlation_supercell)
		pya=phononyaml(filename)
		self.pya=pya
		self.testk=k
		if not exists('NMA'):mkdir('NMA')
		self.specialk(pya,map(float,k))
	def lifenma(self,filename="qpoints/qpoints.yaml",correlation_supercell=[10,10,1]):
		self.pfactor=self.getpfactor(correlation_supercell)
		pya=phononyaml(filename)
		f=open('lifenma.txt','w')
		c=h5py.File('lifenma.h5')
		f.write('kx\tky\tkz\tfreq\tw0\tscatter\ttao\n')
		for iqp in range(pya.nqpoint):
			for ibr in range(pya.nbranch):
				node='/%s/%s'%(iqp,ibr)
				print node
				if not node in c:
					k=pya.qposition(iqp)
					freq=pya.frequency(iqp,ibr)
					vec=pya.atoms(iqp,ibr)
					v=self.calculateLife((k,freq,vec),iqp,ibr)
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
	def fitpart(self,x,q,n=None):
		if n is None:
			n=q.argmax()
		p0 = np.array([x[n], 0.01, q[n]]) #Initial guess
		p=self.fitLife(x,q,p0)
		return p
	def lifeSED(self,k,na,iqp,pya):

		q=self.calculateSED(k,na)
		x=self.freq
		nbranch=3*na
		w=[pya.frequency(iqp,ibr) for ibr in range(nbranch)]
		w1=[]
		df=1.0/(2.0*self.timestep)/(self.totalStep/2)
		span=int(1.0/df)
		for i in range(nbranch):
			node='/%s/%s'%(iqp,i)
			print node
			# get the 1Hz neighbor of the original w,because dw couldn't be too large
			ori=int(w[i]/df)
			low=max(ori-1*span,0)
			hi=min(ori+1*span,len(q)-1)
			filter=range(low,hi)
			qu=q[filter]
			xu=x[filter]
			qs=self.lowess(xu,qu)
			if False:
				from scipy import signal
				peaks, = signal.argrelmax(qs, order=span/3)
				center=peaks[np.abs(xu[peaks]-w[i]).argmin()]
				
			p=self.fitpart(xu,qs)
			w1.append(p)
			#q[filter]=0.0
		w1=np.array(w1)
		
		#o=w1[:,0].argsort()
		#w1=w1[o]
		v=['\t'.join(map(str,list(k)+[w[i],p[0],p[1],1.0/p[1]]) ) for i,p in enumerate(w1)]
		c='\n'.join(v)
		print "[sed]",c
		return c

	def lowess(self,x, y, f=.1):
		f=np.sqrt(10)*self.timestep
		n = len(x)
		r = int(np.ceil(f*n))
		h = [np.sort(np.abs(x - x[i]))[r] for i in range(n)]
		w = np.clip(np.abs((x[:, None] - x[None, :]) / h), 0.0, 1.0)
		w = (1.0 - w**3)**3
		y_smooth = np.zeros(n)
		delta = np.ones(n)
		for i in range(n):
			weights = delta * w[:, i]
			b = np.array([np.sum(weights*y), np.sum(weights*y*x)])
			a = np.array([[np.sum(weights), np.sum(weights*x)],
							[np.sum(weights*x), np.sum(weights*x*x)]])
			beta = np.linalg.solve(a, b)
			y_smooth[i] = beta[0] + beta[1]*x[i]

		return y_smooth


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
		
