from numpy.fft import rfft, irfft,fft,ifft
import numpy as np
from aces.graph import series,plot,imshow
from aces.tools import exit,write,to_txt
from aces.velocityh5 import velocity
from scipy.optimize import leastsq
from aces.dos import plot_smooth,plot_dos,plot_vacf,plot_atomdos
from math import pi,sqrt
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
	def dis(self,pos,i,j,offset):
		x=pos[i,0]-pos[j,0]-offset[0]
		y=pos[i,1]-pos[j,1]-offset[1]
		z=pos[i,2]-pos[j,2]-offset[2]
		return sqrt(x*x+y*y+z*z)
	def get_nei(self,r,dr,pos,cell):
		n=len(pos)
		x=int(r/np.linalg.norm(cell[0]))+1
		y=int(r/np.linalg.norm(cell[1]))+1
		z=int(r/np.linalg.norm(cell[2]))+1
		c=[]
		for i in range(n):
			print i
			for j in range(n):
				
				q=[]
				for ix in range(-x,x+1):
					for iy in range(-y,y+1):
						for iz in range(-z,z+1):
							d=self.dis(pos,i,j,ix*cell[0]+iy*cell[1]+iz*cell[2])
							if d>r:continue
							ir=int(d/dr)
							if ir==0:continue
							q.append(ir)
				c.append([i,j,q])
		return c
	def cal_lc(self,r,dr,m):
		from aces.lammpsdata import lammpsdata
		atoms=lammpsdata().set_src('correlation_structure')
		atoms.set_pbc([m.xp,m.yp,m.zp])
		nw=self.totalStep/2+1
		nr=int(r/dr)+1		
		n=len(atoms)
		dis=atoms.get_all_distances(mic=True)
		pos=atoms.positions
		cell=atoms.cell
		symbols=atoms.get_chemical_symbols()
		masses=self.getMassFromLabel(symbols)
		if not exists('g.npz'):
			g=np.zeros([nw,nr],dtype=np.complex)
			ng=np.zeros(nr,dtype=np.int)
			c=self.get_nei(r,dr,pos,cell)
			for i,j,q in c:
				if j==0:
					print "center atom:%d/%d"%(i,n)
				aij=np.einsum('ij,ij->i',self.fourier_atom(i).conjugate(),self.fourier_atom(j))
				aii=np.einsum('ij,ij->i',self.fourier_atom(i).conjugate(),self.fourier_atom(i))
				ajj=np.einsum('ij,ij->i',self.fourier_atom(j).conjugate(),self.fourier_atom(j))
				x=aij/np.sqrt(aii*ajj)*np.sqrt(masses[i]*masses[j])
				for ir in q:
					g[:,ir]+=x
					ng[ir]+=1
			np.savez('g.npz',g=g,ng=ng)
		npz=np.load('g.npz')
		g=npz['g']
		ng=npz['ng']
		rs=(np.arange(nr))*dr
		rs=rs[:-1]
		def aa(x):
			if x>0:return 1.0/x 
			else: return 0.0
		ivr=[aa(x) for x in ng]#1.0/4.0/np.pi/rs**2
		g=np.einsum('ij,j->ij',g,ivr)[:,:-1]
		g=(g.conjugate()*g).real
		imshow(g,'g.png',extent=[0,1,0,1])
		data=[]
		for i in range(0,self.totalStep/2+1,self.totalStep/20):
			data.append([rs,g[i,:],str(i)])
		series(xlabel='length (A)',
		ylabel='Cross Phonon Energy Spectrum',
		datas=data
		,linewidth=2
		,filename='cpes.png')
		cg=g.cumsum(axis=1)
		imshow(cg,'cg.png',extent=[0,1,0,1])
		cg=np.einsum('ij,i->ij',cg,1.0/cg[:,-1])
		data=[]
		for i in range(0,self.totalStep/2+1,self.totalStep/10):
			data.append([rs,cg[i,:],str(i)])
		series(xlabel='length (A)',
		ylabel='Cross Phonon Energy Spectrum',
		datas=data
		,linewidth=2
		,filename='cges.png')
		lc=((cg<0.80).sum(axis=1))*dr
		x=np.linspace(0,1,self.totalStep/2+1)*1/2.0/self.timestep
		plot([x,'Frequency (THz)'],[lc,'Coherence length(A)'],'coherence.png')
		return lc
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
		node,label=self.inith5('velocity_atom',self.totalStep)
		if not self.db[label][id]:
			print 'prepare %s:%d'%(node,id)

			v=self.velocity.atom(id)
			self.db[node][id]=v
			self.db[label][id]=1
		return self.db[node][id]
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
				sphases.append(np.exp(1j*2.0*pi*(freq+1j*.05)*times))
				vecs.append(vec)
			sphases1=[]
			vecs1=[]
			iqp=1
			for ibr in range(pya.nbranch):
				k=pya.qposition(iqp)
				freq=pya.frequency(iqp,ibr)
				vec=pya.atoms(iqp,ibr)
				natom_unitcell=len(vec)
				phase=self.getPhase(k,natom_unitcell).conjugate()
				vec=np.einsum('ij,i->ij',np.tile(vec,[self.natom/natom_unitcell,1]),phase)
				sphases1.append(np.exp(1j*2.0*pi*freq*times))
				vecs1.append(vec)
			self.sphases=sphases
			self.vecs=vecs
			self.sphases1=sphases1
			self.vecs1=vecs1
			p=np.random.rand(100,pya.nbranch)
			p=np.exp(2j*pi*p)*p
			self.rp=p
			self.testp=True
			self.nbranch=pya.nbranch
		p=self.rp		

		sphases=self.sphases
		vecs=self.vecs
		v=np.zeros([self.totalStep,3],dtype=np.complex)
		for ibr in range(self.nbranch):
			vec=vecs[ibr]
			for i in range(3):				
				v[:,i]+=vec[id,i]*sphases[ibr]*p[iseed,ibr]
		#sphases=self.sphases1
		#vecs=self.vecs1
		#for ibr in range(self.nbranch):
		#	vec=vecs[ibr]
		#	for i in range(3):				
		#		v[:,i]+=vec[id,i]*sphases[ibr]*p[iseed,ibr]
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
	def validatek(self,k):
		k=self.kcell.T.dot(k)
		r=self.scell
		return np.einsum('i,ji',k,r)
	def getPhase(self,k,natom_unitcell):
		p=np.exp(self.pfactor.dot(k))
		#print p
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
		
		phase=self.getPhase(k,natom_unitcell)
		phase1=self.getPhase(-np.array(k),natom_unitcell)
		nd=1
		d=totalStep-nd+1
		phi=np.zeros(d)
		#print phase
		for j in range(natom_unitcell):
			q=np.zeros([totalStep,3],dtype=np.complex)
			#q1=np.zeros([totalStep,3],dtype=np.complex)
			for i in range(j,self.natom,natom_unitcell):
				fv=self.velocity_atom(i)
				q+=np.array(fv)*phase[i]
				#q+=np.array(fv)*phase1[i]
			for k0 in range(nd):
				u0=k0+d
				#if u0>len(q):
				#	u0=len(q)
				u1=u0-d
				q0=fft(q[u1:u0],axis=0)
				#q0=q0[:totalStep/2+1]
				phi+=(q0*q0.conjugate()).real.sum(axis=1)
		phi=phi[:d/2+1]/nd
		s=np.zeros_like(phi)
		df=1.0/(2.0*self.timestep)/(self.totalStep/2)
		span=int(1.0/df)
		phi=self.smo(phi,50)
		if self.totalsed:
			x=np.linspace(0,1,d/2+1)*1/2.0/self.timestep

			series(xlabel='Frequency (THz)',
				ylabel='Phonon Energy Spectrum',
				datas=[(x,phi,"origin")]
				,linewidth=1
				,filename='NMA/sed%s.png'%(str(k)),logy=True)
		return phi
	
	def calculateLife(self,eigen,iqp,ibr):
		totalStep=self.totalStep
		k,freq,vec=eigen
		vec0=vec
		#vec[0,0]=1/np.sqrt(2)
		#vec[1,0]=1/np.sqrt(2)*-1j
		#vec[0,0],vec[1,0]=.5*(1.0+1j),1.0/np.sqrt(2)
		natom_unitcell=len(vec)
		phase=self.getPhase(np.array(k),natom_unitcell)
		vec=np.einsum('ij,i->ij',np.tile(vec0.conjugate(),[self.natom/natom_unitcell,1]),phase)
		phase1=self.getPhase(-np.array(k),natom_unitcell)
		#print np.tile(vec.conjugate(),[self.natom/natom_unitcell,1])[:,0]
		vec1=np.einsum('ij,i->ij',np.tile(vec0,[self.natom/natom_unitcell,1]),phase1)
		#print vec[:,0]
		if not self.test:
			q=np.zeros(totalStep,dtype=np.complex)
			q1=np.zeros(totalStep,dtype=np.complex)
			for i in range(self.natom):
				fv=self.velocity_atom(i)
				q+=np.array(fv).dot(vec[i])
				q1+=np.array(fv).dot(vec1[i])
			#q1=ifft(q)
			q=fft(q)
			q1=fft(q1)
			"""
			result = np.correlate(q, q, mode='full', old_behavior=False)
			keXcorr = result[result.size/2:] / result[result.size/2]

			# keFft: kinetic energy FFT
			keFft = fft(keXcorr[:])
			q=keFft
			"""
			q=q+q1
			q=(q*q.conjugate()).real
			q=q[:totalStep/2+1]
		else:
			self.testp=None
			qc=np.zeros(totalStep/2+1)
			nseed=1
			
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
		span=int(.5/df)
		ori=int(freq/df)
		low=max(ori-1*span,0)
		hi=min(ori+1*span,len(q)-1)
		filter=range(low,hi)
		x=x[filter]
		q=q[filter]/1e10
		#q=self.lowess(x,q)
		q=self.smo(q,3)
		p0 = np.array([x[q.argmax()], 0.01, q[q.argmax()]]) #Initial guess
		p=self.fitLife(x,q,p0)
		p=np.abs(p)
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
		self.scell=atoms.cell
		positions=atoms.positions
		from ase.io import read
		atoms=read('POSCAR')
		self.kcell=atoms.get_reciprocal_cell()
		self.cell=atoms.get_cell()
		return 1j*positions.dot(self.kcell.T)*2*pi
	def findk(self,pya,k0):
		iqp=0
		for iqp0 in range(pya.nqpoint):
			k=pya.qposition(iqp0)
			iqp=iqp0
			if np.allclose(k0,k):
				break
		else:
			raise Exception('no k found!')
		return k,iqp
	def specialk(self,pya,k0=[0,0,0]):
		self.partsed=True
		self.totalsed=True
		k,iqp=self.findk(pya,k0)
		#k1,iqp1=self.findk(pya,-np.array(k0))
		q=self.lifeSED(k,pya.natom,iqp,pya)
		print "check orthorgnal"
		print "k' that k'.dot(k) is not zero :"

		natom_unitcell=pya.natom
		phase=self.getPhase(k,natom_unitcell).conjugate()
		for iqp0 in range(pya.nqpoint):
			k1=pya.qposition(iqp0)
			phase1=self.getPhase(k1,natom_unitcell)
			if np.linalg.norm(np.dot(phase1,phase))>.1:
				print k1,np.linalg.norm(np.dot(phase1,phase))
		print "none zero dot product eigenvector sigma pairs with same k are:"
		p=np.zeros([pya.nbranch]*2)
		for ibr in range(pya.nbranch):
			for ibr1 in range(pya.nbranch):
				vec=pya.atoms(iqp,ibr)
				vec1=pya.atoms(iqp,ibr1)

				#vec=np.einsum('ij,i->ij',np.tile(vec,[self.natom/natom_unitcell,1]),phase)
				#vec1=np.einsum('ij,i->ij',np.tile(vec1,[self.natom/natom_unitcell,1]),phase)
				p[ibr,ibr1]=np.einsum('ij,ij',vec.conjugate(),vec1).real
				if abs(p[ibr,ibr1])>.01:print ibr,ibr1
		print "check Born-Karmen condition:"
		for iqp0 in range(pya.nqpoint):
			k1=pya.qposition(iqp0)
			print self.validatek(k1)
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
		#return
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
			#qs=self.lowess(xu,qu)
			if False:
				from scipy import signal
				peaks, = signal.argrelmax(qs, order=span/3)
				center=peaks[np.abs(xu[peaks]-w[i]).argmin()]
				
			p=self.fitpart(xu,qu)
			w1.append(p)
			#q[filter]=0.0
		w1=np.array(w1)
		
		#o=w1[:,0].argsort()
		#w1=w1[o]
		v=['\t'.join(map(str,list(k)+[w[i],p[0],p[1],1.0/p[1]]) ) for i,p in enumerate(w1)]
		c='\n'.join(v)
		print "[sed]",c
		return c
	def smo(self,y,n):
		s=np.zeros_like(y)
		for i in range(0,n):
			p1=np.roll(y,i)
			p1[:i]=0.0
			s+=p1
		phi=s/float(n)
		return phi
	def lowess(self,x, y, f=.1):
		f=np.sqrt(10)*self.timestep
		f=.4
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
		return np.log(self.lorentz(p,x))-np.log(z)
		
