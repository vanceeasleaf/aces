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
from aces.graph import fig,pl
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

		
	def postT(self):
		from aces.graph import fig,pl
		a=np.loadtxt("BTE.KappaTensorVsT_CONV")
		with fig('T_kappa.png',legend=True):
			ts=a[:,0]
			fil=ts<=800
			k1=a[fil,1]
			k2=a[fil,5]
			k3=a[fil,9]
			ts=a[fil,0]
			pl.plot(ts,k1,lw=2,label="${\kappa_{xx}}$")
			pl.plot(ts,k2,lw=2,label="${\kappa_{yy}}$")
			pl.plot(ts,k3,lw=2,label="${\kappa_{zz}}$")
			pl.xlabel("Tempeature (K)")
			pl.ylabel('Thermal Conductivity (W/mK)')
	def oldT(self):
		from aces.graph import fig,pl
		import matplotlib as mpl 
		mpl.rcParams['axes.color_cycle']=['#e24a33','#2A749A','#988ed5']
		ts=np.arange(200,801,100)
		a=[]
		for dir in ts:
			print dir
			x=np.loadtxt("12844."+str(dir)+"K/BTE.kappa_tensor")
			a.append(x)
		a=np.array(a)
		with fig('T_kappa.png',legend=True):
			k1=a[:,1]
			k2=a[:,5]
			k3=a[:,9]
			pl.plot(ts,k1,lw=3,markersize=30,linestyle='--',markeredgecolor='w',marker=".",label="${\kappa_{xx}}$")
			pl.plot(ts,k2,lw=3,markersize=15,linestyle='--',markeredgecolor='w',marker="v",label="${\kappa_{yy}}$")
			pl.plot(ts,k3,lw=3,markersize=15,linestyle='--',markeredgecolor='w',marker="^",label="${\kappa_{zz}}$")
			pl.xlabel("Tempeature (K)")
			pl.ylabel('Thermal Conductivity (W/mK)')
			pl.xlim([200,800])
	def grtao(self):
		cd('T300K')
		#画格林艾森系数与驰豫时间的关系
		w=np.loadtxt('BTE.w_final')[:,1]
		w=np.abs(w)
		q=np.loadtxt(open('../BTE.qpoints'))
		n=len(q)
		w=w.T.reshape([-1,n])
		w=np.einsum('jk->kj',w)
		w.flags.writeable = True
		omega=np.loadtxt('../BTE.omega')/(2.0*np.pi)
		w[omega<omega.flatten().max()*0.005]=float('nan')
		tao=1.0/w+1e-6
		g=np.loadtxt('../BTE.gruneisen')
		with fig("gruneisen_tao.png"):
			pl.semilogy(g.flatten(),tao.flatten(),ls='.',marker='.',color='r',markersize =10)
			pl.ylabel('Relaxation Time (ps)')
			pl.xlabel('Gruneisen Coeffecient')
			pl.xlim([-10,5])
			pl.ylim([0,1e4])
	def postnew(self):
		cd('T300K')
		try:
			df=pd.read_csv("BTE.kappa_scalar",sep=r"[ \t]+",header=None,names=['step','kappa'],engine='python');
			ks=np.array(df['kappa'])
			plot((np.array(df['step']),'Iteration Step'),(ks,'Thermal Conductivity (W/mK)'),'kappa_scalar.png',grid=True,linewidth=2)
		except Exception as e:
			print e

		try:
			df=pd.read_csv("BTE.cumulative_kappa_scalar",sep=r"[ \t]+",header=None,names=['l','kappa'],engine='python');
			ks=np.array(df['kappa'])
			plot((np.array(df['l']),'Cutoff Mean Free Path for Phonons (Angstrom)'),(ks,'Thermal Conductivity (W/mK)'),'cumulative_kappa_scalar.png',grid=True,linewidth=2,logx=True)	
		except Exception as e:
			print e
		try:
			omega=np.loadtxt('../BTE.omega')/(2.0*np.pi)
			kappa=np.loadtxt('BTE.kappa')[-1,1:]
			kappa=np.einsum('jji',kappa.reshape([3,3,-1]))/3.0
			plot((np.arange(len(omega[0])),'Band'),(kappa,'Thermal Conductivity (W/mK)'),'kappa_band.png',grid=True,linewidth=2)
			plot((np.arange(len(omega[0])),'Band'),(kappa.cumsum(),'Thermal Conductivity (W/mK)'),'cumulative_kappa_band.png',grid=True,linewidth=2)
		except Exception as e:
			print e
		try:
			w=np.loadtxt('BTE.w_final')[:,1]
			w=np.abs(w)
			q=np.loadtxt(open('../BTE.qpoints'))
			n=len(q)
			w=w.T.reshape([-1,n])
			w=np.einsum('jk->kj',w)
			w.flags.writeable = True
			print w.shape,omega.shape
			w[omega<omega.flatten().max()*0.005]=float('nan')
			plot((omega.flatten(),'Frequency (THz)'),(w.flatten(),'Scatter Rate (THz)'),'scatter_freq.png',grid=True,scatter=True,logy=True)
			tao=1.0/w+1e-6
			plot((omega.flatten(),'Frequency (THz)'),(tao.flatten(),'Relaxation Time (ps)'),'tao_freq.png',grid=True,scatter=True,logy=True)
			to_txt(['freq','tao'],np.c_[omega.flatten(),tao.flatten()],'tao_freq.txt')
		except Exception as e:
			print e
		try:
	
			v=np.loadtxt(open('../BTE.v'))
			q=np.loadtxt(open('../BTE.qpoints'))
			n=len(q)
			v=v.T.reshape([3,-1,n])
			v=np.einsum('ijk->kji',v)
			v=np.linalg.norm(v,axis=-1)
			plot((omega.flatten(),'Frequency (THz)'),(v.flatten(),'Group Velocity (nm/ps)'),'v_freq.png',grid=True,scatter=True)
			to_txt(['freq','vg'],np.c_[omega.flatten(),v.flatten()],'v_freq.txt')
	
			
			l=v*tao
			plot((omega.flatten(),'Frequency (THz)'),(l.flatten(),'Mean Free Path (nm)'),'lamda_freq.png',grid=True,scatter=True)
			to_txt(['freq','mfp'],np.c_[omega.flatten(),l.flatten()],'lamda_freq.txt')
		except Exception as e:
			print e
		try:
			g=np.loadtxt('../BTE.gruneisen')
			plot((omega.flatten(),'Frequency (THz)'),(g.flatten(),'Gruneisen'),'gruneisen_freq.png',grid=True,scatter=True)
			with fig('gruneisen_freq.png'):
				pl.scatter(omega.flatten(),g.flatten(),marker='.',color='r',s =50)
				pl.xlabel('Frequency (THz)')
				pl.ylabel('Gruneisen Coeffecient')
				#pl.grid(True)
				pl.xlim([0,omega.max()])
				pl.ylim([-10,5])
				#pl.tick_params(axis='both', which='major', labelsize=14)
			to_txt(['freq','gruneisen'],np.c_[omega.flatten(),g.flatten()],'gruneisen_freq.txt')
			g=np.loadtxt('../BTE.P3')

			with fig('p3_freq.png'):
				pl.scatter(omega.flatten(),g.flatten()*1e6,marker='.',color='r',s =50)
				pl.xlabel('Frequency (THz)')
				pl.ylabel('P3 $(\\times 10^{-6})$')
				#pl.grid(True)
				pl.xlim([0,omega.max()])
				pl.ylim([0,g.max()*1e6])

			to_txt(['freq','p3'],np.c_[omega.flatten(),g.flatten()],'p3_freq.txt')
		except Exception as e:
			print e
		try:	
			q=np.loadtxt(open('../BTE.qpoints'))
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
		cd('..')
	def vtao(self):
		#group velocity vs. tao using old version of shengbte
		w=np.loadtxt('BTE.w_final')
		w=np.abs(w)
		omega=np.loadtxt('BTE.omega')/(2.0*np.pi)
		w[omega<omega.flatten().max()*0.005]=float('nan')
		tao=1.0/w+1e-6
		v=np.loadtxt(open('BTE.v'))
		n,m=v.shape
		v=v.reshape([n,3,m/3])
		v=np.linalg.norm(v,axis=1)
		l=v*tao
		l[l<1e-6]=None
		with fig('tao_v.png'):
			pl.semilogy(v.flatten(),tao.flatten(),linestyle='.',marker='.',color='r',markersize =5)
			pl.xlabel('Group Velocity (nm/ps)')
			pl.ylabel('Relaxation Time (ps)')
			pl.grid(True)
		with fig('tao_l.png'):
			pl.loglog(l.flatten(),tao.flatten(),linestyle='.',marker='.',color='r',markersize =5)
			pl.xlabel('Mean Free Path (nm)')
			pl.ylabel('Relaxation Time (ps)')
			pl.grid(True)
		with fig('v_l.png'):
			pl.semilogy(v.flatten(),l.flatten(),linestyle='.',marker='.',color='r',markersize =5)
			pl.xlabel('Group Velocity (nm/ps)')
			pl.ylabel('Mean Free Path (nm)')
			pl.grid(True)
	def getGrid(self):
		s=shell_exec("grep ngrid CONTROL")
		from scanf import sscanf
		grids=sscanf(s,"ngrid(:)=%d %d %d")
		return grids
	def getQ(self):
		from ase import io
		atoms=io.read('../POSCAR')
		rcell=atoms.get_reciprocal_cell()
		grid=self.getGrid()
		q0=[]
		for ii in range(grid[0]):
			for jj in range(grid[1]):
				for kk in range(grid[2]):
					k=[float(ii)/grid[0]-.5,float(jj)/grid[1]-.5,float(kk)/grid[2]-.5]
					#q0.append(np.einsum('ij,i',rcell,k))
					q0.append(k)
		return np.array(q0)

	def getQFmap(self):
		qpoints_full=np.loadtxt('BTE.qpoints_full')
		qpoints=np.loadtxt('BTE.qpoints')
		ids=qpoints_full[:,1].astype(np.int)
		idx=qpoints[:,0].astype(np.int)
		a={}
		for i,id in enumerate(idx):
			a[id]=i
		u=np.array([a[i] for i in ids])
		return u
	def taoth(self):
		#tao vs. direction in xy plane using old version of shengbte
		w=np.loadtxt('BTE.w_final')
		w=np.abs(w)
		omega=np.loadtxt('BTE.omega')/(2.0*np.pi)
		w[omega<omega.flatten().max()*0.005]=float('nan')
		tao=1.0/w+1e-6
		tao[tao>10000]=0
		tao=np.nan_to_num(tao)
		u=self.getQFmap()
		tao=tao[u]
		#为了限制q点在BZ,必须自己重新来
		#qpoints_full=np.loadtxt('BTE.qpoints_full')
		#q=qpoints_full[:,-3:]
		q=self.getQ()
		with fig('tao_th.png'):
			ax = pl.subplot(111, projection='polar')
			N=100
			th=np.linspace(0,1,N)*np.pi*2.0-np.pi
			r=np.zeros_like(th)
			r1=np.zeros_like(th)
			theta=np.arctan2(q[:,1],q[:,0])
			for i in np.arange(1):
				for j,tt in enumerate(th):
					if j==len(th)-1:						
						fil=(theta>=tt)
					else:
						fil=(theta>=tt) *(theta<th[j+1])
						r[j]=np.nan_to_num(tao[fil].mean())
						r1[j]=np.nan_to_num(fil.sum())
				c = pl.plot(th, r,lw=2)
				#pl.plot(th, r1,lw=2)
				#c.set_alpha(0.75)
				#pl.semilogy(q[:,0].flatten(),tao[:,i].flatten(),linestyle='.',marker='.',color='r',markersize =5)
			pl.grid(True)
	def postsheng(self):
		try:
			df=pd.read_csv("BTE.kappa_scalar",sep=r"[ \t]+",header=None,names=['step','kappa'],engine='python');
			ks=np.array(df['kappa'])
			plot((np.array(df['step']),'Iteration Step'),(ks,'Thermal Conductivity (W/mK)'),'kappa_scalar.png',grid=True,linewidth=2)
		except Exception as e:
			print e	

		try:
			df=pd.read_csv("BTE.cumulative_kappa_scalar",sep=r"[ \t]+",header=None,names=['l','kappa'],engine='python');
			ks=np.array(df['kappa'])
			plot((np.array(df['l']),'Cutoff Mean Free Path for Phonons (Angstrom)'),(ks,'Thermal Conductivity (W/mK)'),'cumulative_kappa_scalar.png',grid=True,linewidth=2,logx=True)	
		except Exception as e:
			print e	
		try:
			omega=np.loadtxt('BTE.omega')/(2.0*np.pi)
			kappa=np.loadtxt('BTE.kappa')[-1,1:]
			kappa=np.einsum('jji',kappa.reshape([3,3,-1]))/3.0
			plot((np.arange(len(omega[0])),'Band'),(kappa,'Thermal Conductivity (W/mK)'),'kappa_band.png',grid=True,linewidth=2)
			plot((np.arange(len(omega[0])),'Band'),(kappa.cumsum(),'Thermal Conductivity (W/mK)'),'cumulative_kappa_band.png',grid=True,linewidth=2)
		except Exception as e:
			print e	
		try:
			w=np.loadtxt('BTE.w_final')
			w=np.abs(w)
			w[omega<omega.flatten().max()*0.005]=float('nan')
			plot((omega.flatten(),'Frequency (THz)'),(w.flatten(),'Scatter Rate (THz)'),'scatter_freq.png',grid=True,scatter=True,logy=True)
			tao=1.0/w+1e-6
			with fig('tao_freq.png'):
				pl.semilogy(omega.flatten(),tao.flatten(),linestyle='.',marker='.',color='r',markersize =5)
				pl.xlabel('Frequency (THz)')
				pl.ylabel('Relaxation Time (ps)')
				pl.grid(True)
				pl.xlim([0,omega.max()])
				#pl.ylim([0,tao.flatten().max()])
			to_txt(['freq','tao'],np.c_[omega.flatten(),tao.flatten()],'tao_freq.txt')
		except Exception as e:
			print e	
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
			print e	
		try:	
			l=v*tao
			l[l<1e-6]=None
			plot((omega.flatten(),'Frequency (THz)'),(l.flatten(),'Mean Free Path (nm)'),'lamda_freq.png',grid=True,scatter=True,logy=True,logx=True,xmin=0)
			to_txt(['freq','mfp'],np.c_[omega.flatten(),l.flatten()],'lamda_freq.txt')
		except Exception as e:
			print e	
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
			print e	
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
	def runold(self):
		#Thermal conductivity calculation
		m=self.m
		print "START SHENGBTE..."
		passthru(config.mpirun+" %s "%(m.nodes*m.procs)+config.sheng)
	def runsheng(self):
		#Thermal conductivity calculation
		m=self.m
		print "START SHENGBTE..."
		passthru(config.mpirun+" %s "%(m.nodes*m.procs)+config.shengbte)
		

	def kkappa(self):
		dirs=ls('shengold*')
		from aces.scanf import sscanf
		from aces.graph import fig,pl
		import matplotlib as mpl 
		mpl.rcParams['axes.color_cycle']=['#e24a33','#2A749A','#988ed5']
		us=[]
		for d in dirs:
			f=shell_exec('grep ngrid %s/CONTROL'%d)
			ks=sscanf(f,"	ngrid(:)=%d %d %d")
			f=np.loadtxt('%s/BTE.kappa_tensor'%d)
			print ks
			if len(f.shape)==2:
				x=f[-1]
			else:
				x=f
			x=x[1:].reshape([3,3])
			print x
			us.append([ks,x])
		
		with fig('kkappa_64nn.png',legend=True):
			p1=filter(lambda u: u[0][0]==64,us)
			k1=[]
			k2=[]
			k3=[]
			ks=[]
			for u in p1:
				ks.append(u[0][1])
				k1.append(u[1][0,0])
				k2.append(u[1][1,1])
				k3.append(u[1][2,2])
			f=np.argsort(ks)
			ks=np.array(ks)[f];k1=np.array(k1)[f];k2=np.array(k2)[f];k3=np.array(k3)[f];
			pl.plot(ks,k1,markersize=30,linestyle='--',markeredgecolor='w',marker=".",label="${\kappa_{xx}}$")
			pl.plot(ks,k2,markersize=15,linestyle='--',markeredgecolor='w',marker="v",label="${\kappa_{yy}}$")
			pl.plot(ks,k3,markersize=15,linestyle='--',markeredgecolor='w',marker="^",label="${\kappa_{zz}}$")
			pl.ylim([0,0.35])
			pl.xlim([0,np.array(ks).max()+1])
			pl.xlabel("$Nq_y$ and $Nq_z$")
			pl.ylabel("Themal Conductivity (W/mK)")
		with fig('kkappa_n44.png',legend=True):
			p1=filter(lambda u: u[0][1]==4,us)
			k1=[]
			k2=[]
			k3=[]
			ks=[]
			for u in p1:
				ks.append(u[0][0])
				k1.append(u[1][0,0])
				k2.append(u[1][1,1])
				k3.append(u[1][2,2])
			f=np.argsort(ks)
			ks=np.array(ks)[f];k1=np.array(k1)[f];k2=np.array(k2)[f];k3=np.array(k3)[f];
			pl.plot(ks,k1,markersize=30,linestyle='--',markeredgecolor='w',marker=".",label="${\kappa_{xx}}$")
			pl.plot(ks,k2,markersize=15,linestyle='--',markeredgecolor='w',marker="v",label="${\kappa_{yy}}$")
			pl.plot(ks,k3,markersize=15,linestyle='--',markeredgecolor='w',marker="^",label="${\kappa_{zz}}$")
			pl.ylim([0,0.35])
			pl.xlim([0,np.array(ks).max()+100])
			pl.xlabel("$Nq_x$")
			pl.ylabel("Themal Conductivity (W/mK)")
	def kmfp(self):
		def ff(p,x):
			#return p[0]*(1.0-np.exp(-x**p[2]/p[1]))
			return 1.0/(p[1]/x+1/p[0])-p[2]
			#return p[0]*p[1]**x
		
		def fit(x,z,p0,tt):
			def errorfunc(p,x,z):
				return tt(p,x)-z
			from scipy.optimize import leastsq
			solp, ier = leastsq(errorfunc, 
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
		dirs=ls('shengold*')
		from aces.scanf import sscanf
		from aces.graph import fig,pl
		us=[]
		for d in dirs:
			f=shell_exec('grep ngrid %s/CONTROL'%d)
			ks=sscanf(f,"	ngrid(:)=%d %d %d")
			if(ks[1]!=4):continue
			f=np.loadtxt('%s/BTE.cumulative_kappa_scalar'%d)
			us.append([ks,f])
		with fig('reduce_mfp.png',legend=True,ncol=1):
			for i,u in enumerate(us):
				if i <3:continue
				ks,f=u
				x,y=f[:,0],f[:,1]
				pl.semilogx(x,y,label="Nx= %d "%ks[0],linewidth=2)
			ks,f=us[-1]
			x,y=f[:,0],f[:,1]
			#fil=(x>0)
			#p=fit(x[fil],y[fil],[1,1,1],ff)
			#y1=ff(p,x)
			#pl.semilogx(x,y1,label="fit of Nx= %d "%ks[0],linewidth=2)
			pl.xlabel('Cutoff Mean Free Path for Phonons (Angstrom)')
			pl.ylabel('Thermal Conductivity (W/mK)')
			pl.grid(True)
		with fig('kappa_inv_mpf_inv.png',legend=True,ncol=1):
			ks,f=us[-1]
			fil=x>.5
			x,y=f[fil,0],f[fil,1]
			xx=1/x;yy=1/y
			pl.plot(xx,yy,linewidth=3,c='red',label="Nx=1024")
			def ll(p,x):
				return p[0]*x+p[1]
			fil=xx>xx.max()/4
			p=fit(xx[fil],yy[fil],[1,1,1],ll)
			pl.plot(xx,ll(p,xx),lw=3,ls='dashed',label="Fitted")
			pl.xlabel('1/L (1/Angstrom)')
			pl.ylabel('$1/\\kappa_L$ (mK/W)')
			pl.grid(True)

