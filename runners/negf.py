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
class runner(Runner):
	def creatmini(self,m):
		print 'creatmini'
		mkdir('minimize')
		cd('minimize')
		minimize_input(m)
		write(time.strftime('%Y-%m-%d %H:%M:%S'),'done')
		cd('..')
	def generate(self):
		fccenter,centeratoms=self.generateCenter()
		fclead,leadatoms=self.generateLead()
		
	def generateCenter(self):
		m=self.m
		mkdir('center')
		cd('center')
		m.supercell=[1,1,1]
		m.phofc=True
		PRunner(m).generate()
		cd('..')
		fccenter=self.readfc('center/FORCE_CONSTANTS')
		centeratoms=read('center/POSCAR')
		fccenter=self.nomalizeFC(fccenter,m,centeratoms)
		return fccenter,centeratoms
	def generateLead(self):
		m=self.m
		mkdir('lead')
		cd('lead')
		s=im('aces.materials.%s'%m.species)
		mm=s.structure(dict(latx=4,laty=1,latz=1))
		mm.supercell=[1,1,1]
		mm.phofc=True
		mm.cores=m.cores
		mm.home=pwd()
		assert mm.home!=''
		self.creatmini(mm)	
		PRunner(mm).generate()	
		cd('..')
		fclead=self.readfc('lead/FORCE_CONSTANTS')
		leadatoms=read('lead/POSCAR')
		fclead=self.nomalizeFC(fclead,mm,leadatoms)
		return fclead,leadatoms
	def rearangefc(self,fc,atoms):

	def readfc(self,filename='FORCE_CONSTANTS'):
		f=open(filename)
		line=f.next()
		natom=int(line)
		fc=np.zeros([natom,natom,3,3])
		for i in range(natom):
			for j in range(natom):
				f.next()
				for k in range(3):
					fc[i,j,k]=map(float,f.next().split())
		return fc
	def nomalizeFC(self,fc,m,atoms):
		natom=len(atoms)
		newfc=np.zeros([natom,natom,3,3])
		masses=m.getMassFromLabel(self.atoms.get_chemical_symbols())
		for i in range(natom):
			for j in range(natom):
				m1=masses[i]
				m2=masses[j]
				newfc[i,j]=-1.0/np.sqrt(m1*m2)*fc[i,j]
		return newfc
