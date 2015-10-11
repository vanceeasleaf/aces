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
from aces.runners.phonopy import runner as Runner
import pandas as pd
class runner(Runner):
			
	def force_constant3(self):
		m=self.m
		cmd='find dirs/dir_3RD.* -name vasprun.xml |sort -n|'+config.thirdorder+" reap %s %f "%(m.dim,m.shengcut/10.0)
		passthru(cmd)

	def generate_supercells3(self):
		m=self.m
		#generate supercells
		cmd=config.thirdorder+"sow %s %f "%(m.dim,m.shengcut/10.0)
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
"""%(len(elements),len(atoms),' '.join(map(str,m.kpoints)))
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
"""%(' '.join(map(str,cell[0])),
			' '.join(map(str,cell[1])),
			' '.join(map(str,cell[2])),
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
	def post(self):
		cd('SHENG')
		df=pd.read_csv("BTE.kappa_scalar",sep=r"[ \t]+",header=None,names=['step','kappa'],engine='python');
		plot((np.array(df['step']),'Iteration Step'),(np.array(df['kappa']),'Thermal Conductivity (W/mK)'),'kappa_scalar.png',grid=True,linewidth=2)
		df=pd.read_csv("BTE.cumulative_kappa_scalar",sep=r"[ \t]+",header=None,names=['l','kappa'],engine='python');
		plot((np.array(df['l']),'Cutoff Mean Free Path for Phonons (Angstrom)'),(np.array(df['kappa']),'Thermal Conductivity (W/mK)'),'cumulative_kappa_scalar.png',grid=True,linewidth=2,logx=True)
		
		omega=np.loadtxt(open('BTE.omega'))/(2.0*np.pi)
		w=np.loadtxt(open('BTE.w'))

		plot((omega.flatten(),'Frequency (THz)'),(w.flatten(),'Scatter Rate (THz)'),'scatter_freq.png',grid=True,scatter=True)
		tao=1.0/w+1e-6
		plot((omega.flatten(),'Frequency (THz)'),(tao.flatten(),'Relaxation Time (ps)'),'tao_freq.png',grid=True,scatter=True,logx=True,logy=True)
		if not exists('relaxtime'):mkdir('relaxtime')
		cd('relaxtime')
		for i,om in enumerate(omega[:6]):
			print "q : ",i
			plot((om,'Frequency (THz)'),(tao[i],'Relaxation Time (ps)'),'tao_freq_q%d.png'%i,grid=True,scatter=True,logx=True,logy=True)
		cd('..')
		v=np.loadtxt(open('BTE.v'))
		n,m=v.shape
		v=v.reshape([n,m/3,3])
		v=np.linalg.norm(v,axis=2)
		plot((omega.flatten(),'Frequency (THz)'),(v.flatten(),'Group Velocity (nm/ps)'),'v_freq.png',grid=True,scatter=True)
		
		l=v*tao
		plot((omega.flatten(),'Frequency (THz)'),(l.flatten(),'Mean Free Path (nm)'),'lamda_freq.png',grid=True,scatter=True)
		
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
		cd('..')
	def generate(self):
		m=self.m
		self.minimizePOSCAR()
		#cp('minimize/POSCAR','.')
		mkdir('secondorder')
		cd('secondorder')
		cp('../POSCAR','.')
		self.generate_supercells()
		files=shell_exec("ls *-*").split('\n')
		assert len(files)>0
		self.getvasprun(files)
		self.force_constant(files)
		cd('..')
		mkdir('thirdorder')
		cd('thirdorder')
		cp('../POSCAR','.')
		self.generate_supercells3()
		files=shell_exec("ls 3RD.*.*|sort -n").split('\n')
		assert len(files)>0
		self.getvasprun(files)
		self.force_constant3()
		cd('..')
		mkdir('SHENG')
		cd('SHENG')
		cp('../secondorder/FORCE_CONSTANTS','FORCE_CONSTANTS_2ND')
		cp('../thirdorder/FORCE_CONSTANTS_3RD','.')
		self.getControl()
		#Thermal conductivity calculation
		print "START SHENGBTE..."
		passthru(config.mpirun+" %s "%m.cores+config.shengbte)
		

