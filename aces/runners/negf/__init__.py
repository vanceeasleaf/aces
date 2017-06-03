#encoding:utf8
from aces.tools import *
import aces.config as config
from ase.io import read
from ase.io.vasp import write_vasp
from aces.runners import Runner
from aces.graph import plot,series
import numpy as np
from aces.runners.minimize import minimize as minimize_input
from aces.runners.phonopy import runner as PRunner
from importlib import import_module as im
import time
from aces.f import readfc2
hbar=6.6260755e-34/3.14159/2.0
kb=1.3806488e-23
def BE(w,T):
	w=np.array(w)
	t= hbar*w/kb/T
	#return np.exp(-t)
	return 1.0/(np.exp(t)-1.0000001)
class runner(Runner):
	def creatmini(self,m):
		print 'creatmini'
		m.home=pwd()
		assert m.home!=''
		mkdir('minimize')
		cd('minimize')
		minimize_input(m)
		write(time.strftime('%Y-%m-%d %H:%M:%S'),'done')
		cd('..')
		return m.dump2POSCAR(m.home+'/minimize/range')
	def test(self):
		dm=.1
		omega=np.arange(dm,60,dm)#THz
		factor=1e12**2*1e-20*1e-3/1.6e-19/6.23e23
		energies=(omega*2.0*np.pi)**2*factor
		#energies=np.arange(0,10,.01)
		h = -np.array((-2, 1,0, 1, -2,1,0,1,-2)).reshape((3,3))
		h1 = -np.array((-2, 1, 1, -2)).reshape((2,2))
		#x=1.0/np.sqrt(2)
		#h1=h=-np.array((-2,x,0,0,x,-1,x,0,0,x,-2,x,0,0,x,-1)).reshape((4,4))
		#energies = np.arange(-3, 3, 0.1)
		calc = TransportCalculator(h=h, h1=h1, energies=energies,dos=True)
		T = calc.get_transmission()
		#print T
		dos=calc.get_dos()*omega
		plot([omega,'Frequency (THz)'],[T,'Transmission'],'test_green_transmission.png')
		plot([omega,'Frequency (THz)'],[dos,'Phonon Density of State'],'test_green_dos.png')
	def collect(self):		
		leadm=self.preLead()
		fclead=self.fc('lead')
		fccenter=self.fc('center')
		#write(np.around(fc[:,:,0,0],3),'orig_forces')
		n=leadm.hatom
		fccenter[:n,-n:]=0
		fccenter[-n:,:n]=0
		#write(np.around(fccenter[:,:,0,0],3),'fccenter')
		fclead=fclead[:2*n][:2*n]
		fccenter=self.reshape(fccenter)
		fclead=self.reshape(fclead)	
		return fccenter,fclead
	def gettrans(self):
		print("Reading in force constants...")
		if not exists("fcbin.npz"):
			fccenter,fclead=self.collect()
			np.savez("fcbin.npz",fccenter=fccenter,fclead=fclead)
			print("Caching force constans")
		import os
		m=self.m
		os.system(config.mpirun+" "+str(m.cores)+" ae trans_cal >trans_cal.out")
	def reduce(self):
		files=ls("tmp/result.txt*")
		omega=[]
		trans=[]
		dos=[]
		for file in files:
			result=np.loadtxt(file,skiprows=1)
			omega.append(result[:,0])
			trans.append(result[:,1])
			dos.append(result[:,2])

		omega=np.array(omega).flatten().T
		f=omega.argsort()
		omega=omega[f]
		trans=np.array(trans).flatten().T[f]
		dos=np.array(dos).flatten().T[f]
		to_txt(['omega','trans','dos'],np.c_[omega,trans,dos],'result.txt')
	def generate(self):
		self.m.xp=1
		leadm=self.preLead()
		self.phonopy('lead',leadm)
		centerm=self.preCenter()
		self.phonopy('center',centerm)
		self.gettrans()
		self.post()
	def post(self):
		#1eV = 8049 cm^(-1) => 1000emV=8049 cm-1 => cm-1/meV=1000/8049
		#1cm^(-1) = 3 * 10^(10) hz =>Hz*cm=1/3e10
		#a cm^-1=b THz =>a=b *1e12 Hz*cm
		#a meV = b cm^-1 => a = b cm-1/meV
		#omcm=omega*521.471？
		result=np.loadtxt("result.txt",skiprows=1)
		omega=result[:,0]
		trans =result[:,1]
		dos   =result[:,2]
		omcm=omega*1e12*1/3e10
		omme=omcm *1e12*6.6260755e-34/1.6e-19*1000
		w=omega*1e12*2.0*np.pi
		T=self.m.T
		centerm=self.preCenter()
		V=np.linalg.det(centerm.atoms.cell)
		c=hbar*w*(BE(w,T+0.005)-BE(w,T-0.005))*100.0/V*1e30
		j=c*trans/2.0/np.pi
		dm=omega[1]-omega[0]
		kappa=j.cumsum()*dm
		to_txt(['Frequency (THz)',
				'Frequency (cm^-1)',
				'Frequency (meV)',
				'Phonon Transmission',
				'Phonon Density of State',
				'Mode Capacity (J/m^3/K)',
				'Mode Thermal Conductance (W/m^2/K)',
				'Accumulate Thermal Conductance (W/m^2/K)'
			],
		np.c_[omega,omcm,omme,trans,dos,c,j,kappa],'transmission.txt')
		
		f=np.loadtxt('transmission.txt',skiprows=1)
		plot([f[:,0],'Frequency (THz)'],[f[:,4],'Phonon Density of State'],'green_dos.png')
		plot([f[:,0],'Frequency (THz)'],[f[:,3],'Phonon Transmission'],'green_transmission.png')
		plot([f[:,0],'Frequency (THz)'],[f[:,6],'Mode Thermal Conductance (W/m^2/K)'],'green_mode_conductance.png')
	def reshape(self,fc):
		n,m=fc.shape[:2]
		fc=np.einsum('ikjl',fc).reshape([n*3,m*3])
		return fc
	def fc(self,dir):
		fc=readfc2(dir+'/FORCE_CONSTANTS')
		atoms=read(dir+'/POSCAR')
		fc=self.nomalizeFC(fc,atoms)		
		
		fc=self.rearangefc(fc,atoms,dir)

		return fc
	def phonopy(self,dir,mm):
		if exists(dir+'/FORCE_CONSTANTS'):
			return
		mkcd(dir)	
		self.creatmini(mm)
		PRunner(mm).generate()
		cd('..')
		
	def preCenter(self):
		m=self.m
		import device.device as s
		leadm=self.preLead()
		u=s.Device(m,leadm,leadm)
		u.cores=m.cores
		u.__dict__=dict(m.__dict__,**u.__dict__)
		return u

	def preLead(self):
		m=self.m
		s=im('aces.materials.%s'%m.leads)
		lat=m.leadlat
		mm=s.structure(dict(latx=lat[0],laty=lat[1],latz=lat[2],xp=1,yp=1,zp=1))
		mm.dimension=m.dimension
		import device.lead as s
		u=s.Lead(mm)
		#u.cores=m.cores
		u.__dict__=dict(m.__dict__,**u.__dict__)
		return u
	def rearangefc(self,fc,atoms,dir):
		#from aces.f import mapatoms,writefc2
		#pos,order=mapatoms(atoms,old)
		#old.write('old.xyz')
		#atoms[order].write('new.xyz')
		order=np.loadtxt(dir+'/POSCARswap').astype(np.int)
		#writefc2(fc[order][:,order],'fc')
		return fc[order][:,order]

	def nomalizeFC(self,fc,atoms):
		from aces.materials import getMassFromLabel
		natom=len(atoms)
		newfc=np.zeros([natom,natom,3,3])
		masses=getMassFromLabel(atoms.get_chemical_symbols())
		for i in range(natom):
			for j in range(natom):
				m1=masses[i]
				m2=masses[j]
				newfc[i,j]=1.0/np.sqrt(m1*m2)*fc[i,j]
		return newfc
