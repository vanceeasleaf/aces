#encoding:utf8
import sys
from devices import nvtDevice,mpDevice,ijDevice,gkDevice
from aces.tools import *
import aces.config as config
from ase.io import read
from ase.io.vasp import write_vasp
from aces.runners import Runner
import numpy as np
from aces.graph import plot,series
from numpy.fft import fft,rfft,irfft
def acf(a):
	length = len(a)
	print "start fourier transform of correlation"
	a = rfft(a,axis=0)
	print "start inverse fourier transform of correlation"
	c = irfft(a*a.conjugate(),axis=0)/length
	return c
class Hook:
	def __init__(self):
		self.labels={}
	def addAction(self,label,function):
		if not self.labels.has_key(label):
			self.labels[label]=[]
		self.labels[label].append(function)
	def doAction(self,label):
		if not self.labels.has_key(label):return
		for key in self.labels[label]:
			key()
class runner(Runner):
	def generate(self):
		m=self.m
		hook=Hook()
		if m.method=='nvt':
			device=nvtDevice(hook,m)
		elif m.method=='muller':
			device=mpDevice(hook,m)
		elif m.method=='inject':
			device=ijDevice(hook,m)
		elif m.method=='greenkubo':
			device=gkDevice(hook,m)
		debug("tcfactor="+str(m.tcfactor))
		#settings
		print "units %s"%m.units
		print "dimension 3"
		pbcx=pbcy=pbcz='s'
		if m.xp==1:pbcx='p'
		if m.yp==1:pbcy='p'
		if m.zp==1:pbcz='p'
		print "boundary %s %s %s"%(pbcx,pbcy,pbcz)
		print "atom_style atomic"
		print "read_restart   minimize/restart.minimize"
		print "change_box	all	boundary %s %s %s"%(pbcx,pbcy,pbcz)
		print "lattice fcc 5" #needed to define the regions
		print "thermo %d"%m.dumpRate
		print "thermo_modify     lost warn"
		print "timestep %f"%m.timestep
		#regions and groups
		hook.doAction('region')
		#computes
		print "compute           ke  all  ke/atom"
		print "compute           pe  all  pe/atom"
		print "compute         stress all stress/atom NULL virial"
		print "compute jflux all heat/flux ke pe stress"
		hook.doAction('compute')
		#variables
		hook.doAction('variable')
		#init atoms to T
		print m.masses
		print m.potential
		print "reset_timestep 0"
		print "velocity all create %f %d mom yes rot yes dist gaussian"%(m.T,m.seed)
		if m.dimension==1:
			print "velocity  all set NULL 0.0 0.0 units box"
		elif m.dimension==2:
			print "velocity  all set NULL NULL 0.0 units box"
		hook.doAction('equ')
		
		#/* 定时输出dump文件并按id排序*/
		if(m.dumpxyz):
			print "dump dump1 all atom %d dump.lammpstrj"%(m.dumpRate)
			print "dump_modify  dump1 sort id"
		print "run %d"%m.equTime
		print "unfix getEqu"
		print "reset_timestep 0"
		hook.doAction('elimination')
		print "fix    flux_out  all  ave/time  1  %d  %d  c_jflux[1]  c_jflux[2] c_jflux[3] file  flux.txt "%(m.aveRate,m.aveRate)
		hook.doAction('temp')
		hook.doAction('flux')
		hook.doAction('swap')




		#/* 定时输出速度文件用于计算速度关联函数*/
		if(m.dumpv):
			print "dump dump2 all custom %d dump.velocity type vx vy vz"%(m.dumpRate)
			print "dump_modify  dump2 sort id"
		if m.dimension==1:
			print "fix   1d1 all setforce NULL 0. 0."
		print "run	%d"%(m.runTime)

	def runcmd(self):
		return config.mpirun+"  %s "%self.m.cores+config.lammps+" <input  >log.out"
	def post(self):
		m=self.m
		if m.method=="greenkubo":
			if m.fourierTc:
				self.correlation()
			elif m.computeTc:
				self.reduce(1,m.runTime)
			
	def correlation(self):
		m=self.m
		if not exists("jin.npy"):
			print "loading jin.txt"
			j=np.loadtxt('jin.txt',skiprows=2)[:,1]
			np.save('jin.npy',j)
		j=np.load('jin.npy')
		print "loaded"
		n=len(j)
		if not exists("jj0.npy"):
			
			jj0=acf(j)[:n/2+1]
			np.save('jj0.npy',jj0)
		jj0=np.load('jj0.npy')[:m.aveRate]

		xlo,xhi,ylo,yhi,zlo,zhi,lx,ly,lz=m.box
		v=lx*ly*lz
		factor=m.corRate*m.timestep/(v*m.kb*m.T*m.T)*m.zfactor*m.tcfactor
		x=range(len(jj0))
		plot([x,'Correlation Time (ps)'],[jj0,'Heat Flut Correlation Function'],'acf.png')
		plot([x,'Correlation Time (ps)'],[jj0.cumsum()*factor,'Thermal Conductivity (W/mK)'],'kappa.png')
	def profile(self):
		m=self.m
		cd('minimize')
		m.postMini()
		cd('..')
		import profile
		profile.run(**m.__dict__)
	def q(self):
		self.profile()
		self.cal()
	def cal(self):
		import query as qu
		r={}
		r['kappa']=qu.kappa()
		r['nAtom']=qu.nAtom()
		qu.drawStructure()
		r['ineq']=qu.ineq(self.m)
		print(r)
	def reduce(self,n=1,name0=500000):
		name='ac'+str(name0)+'.dat'
		m=self.m
		xs=[]
		for i in range(n):
			dir='../%d/'%i
			print dir
			x=np.loadtxt(dir+name)
			xs.append(x[:,1])
		ts=x[:,0]
		xs=np.array(xs)
		avexs=xs.mean(axis=0)
		np.savetxt('aveac.txt',np.c_[ts,avexs])
		
		plot([ts,'Correlation Time (ps)'],[avexs,'Heat Flut Correlation Function'],'aveac.png')
		datas=[]
		for i in range(n):
			datas.append([ts,xs[i],'b'])
		series(xlabel='Correlation Time (ps)',ylabel='Heat Flut Correlation Function',
			datas=datas,filename='ac.png',legend=False)
		name='tc'+str(name0)+'.dat'
		xs=[]
		for i in range(n):
			dir='../%d/'%i
			print dir
			x=np.loadtxt(dir+name)
			xs.append(x[:,1])
		ts=x[:,0]
		xs=np.array(xs)
		avexs=xs.mean(axis=0)
		np.savetxt('avetc.txt',np.c_[ts,avexs])
		plot([ts,'Correlation Time (ps)'],[avexs,'Thermal Conductivity (W/mK)'],'avetc.png')
		datas=[]
		for i in range(n):
			datas.append([ts,xs[i],'b'])
		series(xlabel='Correlation Time (ps)',ylabel='Thermal Conductivity (W/mK)',
			datas=datas,filename='tc.png',legend=False)

	def fftac(self):
		x=np.loadtxt('aveac.txt')
		N=len(x)/100
		y=np.fft.rfft(x[:,1])[:N/2+1]
		plot([np.arange(N/2+1),'Frequency (Thz)'],[(y*y.conjugate()).real,'Power Spectrum '],'powerspect.png')



