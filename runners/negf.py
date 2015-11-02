#encoding:utf8
from aces.tools import *
import aces.config as config
from ase.io import read
from ase.io.vasp import write_vasp
from aces.binary import pr
from aces.runners import Runner
from aces.graph import plot,series
import numpy as np
from aces.runners.minimize import minimize as minimize_input
from aces.runners.phonopy import runner as PRunner
from importlib import import_module as im
import time
from ase.transport.calculators import TransportCalculator
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
	def generate(self):
		self.m.xp=1
		self.leadm=self.preLead()
		self.phonopy('lead',self.leadm)
		self.centerm=self.preCenter()
		self.phonopy('center',self.centerm)
		
		fclead=self.fc('lead',self.leadm)
		fccenter=self.fc('center',self.centerm)
		dm=.5
		omega=np.arange(dm,60,dm)#THz
		factor=1e12**2*1e-20*1e-3/1.6e-19/6.23e23
		energies=(omega*2.0*np.pi)**2*factor
		tcalc =TransportCalculator(h=fccenter,h1=fclead,h2=fclead,energies=energies,logfile='negf.log',dos=True)
		
		
		
		print 'Calculate Transmission'
		trans=tcalc.get_transmission()
		print 'Calculate Dos'
		dos=tcalc.get_dos()*omega
		print 'Calculate Thermal Conductance'
		#1eV = 8049 cm^(-1) => 1000emV=8049 cm-1 => cm-1/meV=1000/8049
		#1cm^(-1) = 3 * 10^(10) hz =>Hz*cm=1/3e10
		#a cm^-1=b THz =>a=b *1e12 Hz*cm
		#a meV = b cm^-1 => a = b cm-1/meV
		omcm=omega*1e12*1/3e10
		omme=omcm *1e12*6.6260755e-34/1.6e-19*1000
		w=omega*1e12*2.0*np.pi
		T=self.m.T
		V=np.linalg.det(self.centerm.atoms.cell)
		c=hbar*w*(BE(w,T+0.005)-BE(w,T-0.005))*100.0/V*1e30
		j=c*trans/2.0/np.pi
		kappa=j.cumsum()*dm
		to_txt(['Frequency (THz)','Frequency (cm^-1)','Frequency (meV)','Phonon Transmission','Phonon Density of State','Mode Capacity (J/m^3/K)','Mode Thermal Conductance (W/m^2/K)','Accumulate Thermal Conductance (W/m^2/K)'],np.c_[omega,omcm,omme,trans,dos,c,j,kappa],'transmission.txt')
	def fc(self,dir,mm):
		fc=readfc2(dir+'/FORCE_CONSTANTS')
		atoms=read(dir+'/POSCAR')
		fc=self.nomalizeFC(fc,mm,atoms)		
		
		fc=self.rearangefc(fc,atoms,dir)
		n,m=fc.shape[:2]
		fc=np.einsum('ikjl',fc).reshape([n*3,m*3])
		print fc.shape
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
		s=im('aces.device')
		mm=s.Device(m,self.leadm,self.leadm)
		mm.cores=m.cores
		return mm

	def preLead(self):
		m=self.m
		s=im('aces.materials.%s'%m.leads)
		lat=m.leadlat
		mm=s.structure(dict(latx=lat[0],laty=lat[1],latz=lat[2],xp=1,yp=1,zp=1))
		s=im('aces.lead')
		u=s.Lead(mm)
		u.cores=m.cores
		return u
	def rearangefc(self,fc,atoms,dir):
		#from aces.f import mapatoms,writefc2
		#pos,order=mapatoms(atoms,old)
		#old.write('old.xyz')
		#atoms[order].write('new.xyz')
		order=np.loadtxt(dir+'/POSCARswap').astype(np.int)
		#writefc2(fc[order][:,order],'fc')
		return fc[order][:,order]

	def nomalizeFC(self,fc,m,atoms):
		natom=len(atoms)
		newfc=np.zeros([natom,natom,3,3])
		masses=m.getMassFromLabel(atoms.get_chemical_symbols())
		for i in range(natom):
			for j in range(natom):
				m1=masses[i]
				m2=masses[j]
				newfc[i,j]=1.0/np.sqrt(m1*m2)*fc[i,j]
		return newfc
