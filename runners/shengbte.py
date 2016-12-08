#encoding:utf8
from aces.tools import *
import aces.config as config
from ase.io import read
from ase.io.vasp import write_vasp
from aces.binary import pr
from aces.runners import Runner
from aces.UnitCell.unitcell import UnitCell
from aces.graph import plot,series
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pl
import numpy as np
from aces.f import toString
from aces.runners.phonopy import runner as Runner
import pandas as pd
class runner(Runner):
	def fc3(self):
		self.force_constant3();		
	def force_constant3(self):
		cmd='find dirs/dir_3RD.* -name vasprun.xml |sort -n|'+config.thirdorder+" reap"+self.getcut()
		passthru(cmd)
	def getcut(self):
		m=self.m
		cut=str(m.shengcut/10.0)
		if m.shengcut<0:
			cut=str(m.shengcut)
		return " %s %s "%(m.toString(m.supercell3),cut)
	def generate_supercells3(self):
		#generate supercells
		cmd=config.thirdorder+"sow"+self.getcut()
		print cmd
		passthru(cmd)
	
	def getControl(self):
		
		m=self.m
		
		
		f=open('CONTROL','w')
		atoms=read('../POSCAR')#m.atoms
		elements=m.elements
		#shengbte needs nelements <=natoms
		if len(elements)>len(atoms):elements=elements[:len(atoms)]
		allocations="""&allocations
	nelements=%d
	natoms=%d
	ngrid(:)=%s
&end
"""%(len(elements),len(atoms),toString(m.kpoints))
		cell=atoms.cell
		types=m.toString([m.elements.index(x)+1 for x in atoms.get_chemical_symbols()])
		pos=""
		for i,atom in enumerate(atoms):
			pos+="	positions(:,%d)=%s\n"%(i+1,m.toString(atoms.get_scaled_positions()[i]))
		crystal="""&crystal
	lfactor=0.1,
	lattvec(:,1)=%s
	lattvec(:,2)=%s
	lattvec(:,3)=%s
	elements=%s
	types=%s
%s
	scell(:)=%s
&end
"""%(toString(cell[0]),
			toString(cell[1]),
			toString(cell[2]),
			' '.join(map(lambda x: '"'+x+'"',elements)),
			types,
			pos,
			m.dim)
		parameters="""&parameters
	T=%f
	scalebroad=1.0
&end
"""%(m.T)
		flags="""
&flags
	nonanalytic=.TRUE.
	nanowires=.FALSE.
&end  
		"""
		f.write(allocations)
		f.write(crystal)
		f.write(parameters)
		f.write(flags)
		f.close()
	def sca(self,th=0.0):

		qpoints_full=np.loadtxt('BTE.qpoints_full')
		ks=qpoints_full[:,2:4]
		f=self.direction(ks,th)
		ids=qpoints_full[:,1].astype(np.int)[f]
		qpoints=np.loadtxt('BTE.qpoints')
		idx=qpoints[:,0].astype(np.int)
		u=[list(idx).index(i) for i in ids]

		w=np.loadtxt('BTE.w_final')
		omega=np.loadtxt('BTE.omega')/(2.0*np.pi)
		#w[omega<omega.flatten().max()*0.005]=float('nan')
		tao=1.0/w+1e-6

		rt=tao[u,:3]
		rom=omega[u,:3]
		data=[]
		n,m=rom.shape
		for i in range(m):

			data.append([rom[:,i],rt[:,i],'b'])
		series(xlabel='Frequency (THz)',
		ylabel='Relaxation Time (ps)',
		datas=data,
		filename='scaling-%f.png'%th,scatter=True,legend=False,logx=True,logy=True)
	def sca1(self):
		qpoints_full=np.loadtxt('BTE.qpoints_full')
		ks=qpoints_full[:,2:4]
		f=self.norm(ks,2.3)
		ids=qpoints_full[:,1].astype(np.int)[f]
		qpoints=np.loadtxt('BTE.qpoints')
		idx=qpoints[:,0].astype(np.int)
		u=[list(idx).index(i) for i in ids]
		w=np.loadtxt('BTE.w_final')
		omega=np.loadtxt('BTE.omega')/(2.0*np.pi)
		w[omega<omega.flatten().max()*0.005]=float('nan')
		tao=1.0/w+1e-6

		rt=tao[u,:3]
		rom=omega[u,:3]
		data=[]
		n,m=rom.shape
		for i in range(m):

			data.append([rom[:,i],rt[:,i],'b'])
		series(xlabel='Frequency (THz)',
		ylabel='Relaxation Time (ps)',
		datas=data,
		filename='norm.png',scatter=True,legend=False,logx=True,logy=True)
	def norm(self,ks,r):
		filter=np.abs(np.linalg.norm(ks,axis=1)-r)<1
		return filter
	def sca3(self):
		qpoints_full=np.loadtxt('BTE.qpoints_full')
		ks=qpoints_full[:,2:4]
		f=self.kx(ks,2.3)
		ids=qpoints_full[:,1].astype(np.int)[f]
		qpoints=np.loadtxt('BTE.qpoints')
		idx=qpoints[:,0].astype(np.int)
		u=[list(idx).index(i) for i in ids]
		w=np.loadtxt('BTE.w_final')
		omega=np.loadtxt('BTE.omega')/(2.0*np.pi)
		w[omega<omega.flatten().max()*0.005]=float('nan')
		tao=1.0/w+1e-6

		rt=tao[u,:3]
		rom=omega[u,:3]
		data=[]
		n,m=rom.shape
		for i in range(m):

			data.append([rom[:,i],rt[:,i],'b'])
		series(xlabel='Frequency (THz)',
		ylabel='Relaxation Time (ps)',
		datas=data,
		filename='kx.png',scatter=True,legend=False,logx=True,logy=True)
	def kx(self,ks,r):
		filter=np.abs(ks[:,0]-r)<0.25
		return filter
	def sca2(self):
		w=np.loadtxt('BTE.w_final')
		omega=np.loadtxt('BTE.omega')/(2.0*np.pi)
		w[omega<omega.flatten().max()*0.005]=float('nan')
		tao=1.0/w+1e-6
		rt=tao[50:55,:]
		rom=omega[50:55,:]
		data=[]

		n,m=rom.shape
		for i in range(n):

			data.append([rom[i,:],rt[i,:],'b'])
		series(xlabel='Frequency (THz)',
		ylabel='Relaxation Time (ps)',
		datas=data,
		filename='k.png',scatter=True,legend=False,logx=True,logy=True)
	def direction(self,ks,t):
		#find the k points that are in the t direction
		if t<0:
			t=-t
		t=t*np.pi/180.0
		u=np.arctan2(ks[:,1],ks[:,0])
		b=u-t
		b[b>np.pi]-=2.0*np.pi
		b[b<-np.pi]+=2.0*np.pi		
		filter=np.abs(b)<.5*np.pi/180.0
		return filter
	def postsheng(self):
		try:
			df=pd.read_csv("BTE.kappa_scalar",sep=r"[ \t]+",header=None,names=['step','kappa'],engine='python');
			ks=np.array(df['kappa'])
			plot((np.array(df['step']),'Iteration Step'),(ks,'Thermal Conductivity (W/mK)'),'kappa_scalar.png',grid=True,linewidth=2)
		except Exception as e:
			pass

		try:
			df=pd.read_csv("BTE.cumulative_kappa_scalar",sep=r"[ \t]+",header=None,names=['l','kappa'],engine='python');
			ks=np.array(df['kappa'])
			plot((np.array(df['l']),'Cutoff Mean Free Path for Phonons (Angstrom)'),(ks,'Thermal Conductivity (W/mK)'),'cumulative_kappa_scalar.png',grid=True,linewidth=2,logx=True)	
		except Exception as e:
			pass
		try:
			omega=np.loadtxt('BTE.omega')/(2.0*np.pi)
			kappa=np.loadtxt('BTE.kappa')[-1,1:]
			kappa=np.einsum('jji',kappa.reshape([3,3,-1]))/3.0
			plot((np.arange(len(omega[0])),'Band'),(kappa,'Thermal Conductivity (W/mK)'),'kappa_band.png',grid=True,linewidth=2)
			plot((np.arange(len(omega[0])),'Band'),(kappa.cumsum(),'Thermal Conductivity (W/mK)'),'cumulative_kappa_band.png',grid=True,linewidth=2)
		except Exception as e:
			pass
		try:
			w=np.loadtxt('BTE.w_final')
			w=np.abs(w)
			w[omega<omega.flatten().max()*0.005]=float('nan')
			plot((omega.flatten(),'Frequency (THz)'),(w.flatten(),'Scatter Rate (THz)'),'scatter_freq.png',grid=True,scatter=True,logy=True)
			tao=1.0/w+1e-6
			plot((omega.flatten(),'Frequency (THz)'),(tao.flatten(),'Relaxation Time (ps)'),'tao_freq.png',grid=True,scatter=True,logy=True)
			to_txt(['freq','tao'],np.c_[omega.flatten(),tao.flatten()],'tao_freq.txt')
		except Exception as e:
			pass
		"""
		if not exists('relaxtime'):mkdir('relaxtime')
		cd('relaxtime')
		for i,om in enumerate(omega[:6]):
			print "q : ",i
			plot((om,'Frequency (THz)'),(tao[i],'Relaxation Time (ps)'),'tao_freq_q%d.png'%i,grid=True,scatter=True,logx=True,logy=True)
		cd('..')
		"""

		try:
			v=np.loadtxt(open('BTE.v'))
			n,m=v.shape
			v=v.reshape([n,3,m/3])
			v=np.linalg.norm(v,axis=1)
			plot((omega.flatten(),'Frequency (THz)'),(v.flatten(),'Group Velocity (nm/ps)'),'v_freq.png',grid=True,scatter=True)
			to_txt(['freq','vg'],np.c_[omega.flatten(),v.flatten()],'v_freq.txt')
		except Exception as e:
			pass	
		try:	
			l=v*tao
			plot((omega.flatten(),'Frequency (THz)'),(l.flatten(),'Mean Free Path (nm)'),'lamda_freq.png',grid=True,scatter=True)
			to_txt(['freq','mfp'],np.c_[omega.flatten(),l.flatten()],'lamda_freq.txt')
		except Exception as e:
			pass	
		try:	
			q=np.loadtxt(open('BTE.qpoints'))
			qnorm=np.linalg.norm(q[:,-3:],axis=1)
			data=[]
			n,m=w.shape
			for i in range(m):
				data.append([qnorm,w[:,i],'b'])
			series(xlabel='|q| (1/nm)',
			ylabel='Scatter Rate (THz)',
			datas=data,
			filename='branchscatter.png',scatter=True,legend=False,logx=True,logy=True)
		except Exception as e:
			pass	
	def third(self):
		mkdir('thirdorder')
		cd('thirdorder')
		cp('../POSCAR','.')
		self.generate_supercells3()
		
	def vasprun3(self):
		files=shell_exec("ls 3RD.*.*|sort -n").split('\n')
		assert len(files)>0
		self.getvasprun(files)
	def pSecond(self):
		cp('../POSCAR','.')
		self.generate_supercells()
		files=shell_exec("ls *-*").split('\n')
		assert len(files)>0
		self.getvasprun(files)
	def generate(self):
		m=self.m
		self.minimizePOSCAR()
		#cp('minimize/POSCAR','.')
		mkdir('secondorder')
		cd('secondorder')
		self.pSecond()
		self.fc2()
		cd('..')
		self.third()
		self.vasprun3()
		self.force_constant3()
		cd('..')
		self.pSheng()
		self.runsheng()
	def pSheng(self):
		mkdir('SHENG')
		cd('SHENG')
		cp('../secondorder/FORCE_CONSTANTS','FORCE_CONSTANTS_2ND')
		cp('../thirdorder/FORCE_CONSTANTS_3RD','.')
		self.getControl()
	def runsheng(self):
		#Thermal conductivity calculation
		m=self.m
		print "START SHENGBTE..."
		passthru(config.mpirun+" %s "%(m.nodes*m.procs)+config.shengbte)
		

