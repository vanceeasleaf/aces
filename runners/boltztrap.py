#encoding:utf8
from aces.runners import Runner
from ase.io.vasp import write_vasp
from ase import io
from aces.tools import *
import numpy as np
from  aces import config
from aces.runners.phonopy import runner as Runner
from aces.f import read_forces,matrixFormat
import vasp2boltz
from ase.lattice.spacegroup import Spacegroup
from aces.graph import fig,pl
class runner(Runner):
	def generate(self):
		cp('minimize/CONTCAR','POSCAR')
		self.prepareBand()
		self.prepareet()
		self.runet()
	def prepareBand(self):
		mkcd('band')
		#run vasp
		cp('../POSCAR','.')
		self.getVaspRun_vasp()
		#cp('/home1/xggong/zhouy/tcscripts/sis/sis1.1/0/secondorder/dirs/dir_POSCAR-002/EIGENVAL','.')
		#cp('/home1/xggong/zhouy/tcscripts/sis/sis1.1/0/secondorder/dirs/dir_POSCAR-002/OUTCAR','.')
		cd('..')
	def prepareet(self):
		atoms = io.read('POSCAR')
		from pyspglib import spglib
		s=spglib.get_spacegroup(atoms)
		from aces.scanf import sscanf
		sgN=sscanf(s,'%s (%d)')[1]
		sg = Spacegroup(sgN)
		atoms.info = {'spacegroup': sg}
		nelect=self.getNelect(outcar='band/OUTCAR')
		# The remaining part takes care of writing the files required by BoltzTraP.
		bandstructure= vasp2boltz.get_vasp_bandstructure(pathname='band/')
		mkcd('et')
		vasp2boltz.write_bandstructure_boltztrap(bandstructure,filename='et.energy')
		vasp2boltz.write_structure_boltztrap(atoms,filename='et.struct')
		vasp2boltz.write_intrans_boltztrap(n_electrons = nelect,filename='et.intrans')
		s=read('et.energy')
		s=s.replace('HTE','et')
		write(s,'et.energy')
		s=read('et.struct')
		s=s.replace('HTE','et')
		write(s,'et.struct')
		cd('..')
	def runet(self):
		cd('et')
		m=self.m
		passthru(config.x_trans+" BoltzTraP")
		cd('..')
	def showCmd(self):
		print config.x_trans+" BoltzTraP"
	def get_outfermi(self):
		file=ls("*.outputtrans")[0]
		a=shell_exec("grep FermiE %s|tail -1"%file);
		from aces.scanf import sscanf
		a= sscanf(a,"FermiE:  %f.")[0]
		print "E-fermi="+str(a)+"Ry"
		return a
		
	def post(self):
		file=ls("*.trace")[0]
		d=np.loadtxt(file,skiprows=1)
		d[:,0]-=self.get_outfermi()
		d[:,0]*=13.6
		T=np.unique(d[:,1])
		T=[300]
		with fig("Seebeck.png"):
			pl.xlabel("Energy (eV)")
			pl.ylabel("Seebeck Coefficient ($\mu$V/K)")
			for t in [T[-1]]:
				idx=d[:,1]==t
				pl.plot(d[idx,0],d[idx,4],label="T="+str(t))
		with fig("kappa.png"):
			pl.xlabel("Energy (eV)")
			pl.ylabel("Electronic Thermal Conductivity (W/mK)")
			for t in [T[-1]]:
				idx=d[:,1]==t
				pl.plot(d[idx,0],d[idx,7],label="T="+str(t))
		with fig("powerfactor.png"):
			pl.xlabel("Energy (eV)")
			pl.ylabel("$S^{2}\sigma/\\tau $")
			for t in [T[-1]]:
				idx=d[:,1]==t
				S=d[idx,4]
				sigma=d[idx,5]
				pl.plot(d[idx,0],S*S*sigma,label="T="+str(t))
		with fig("sigma.png"):
			pl.xlabel("Energy (eV)")
			pl.ylabel("$sigma/\\tau $")
			for t in [T[-1]]:
				idx=d[:,1]==t
				sigma=d[idx,5]
				pl.plot(d[idx,0],sigma,label="T="+str(t))
		with fig("Rh-n.png"):
			pl.xlabel("Energy (eV)")
			pl.ylabel("$sigma/\\tau $")
			for t in [T[-1]]:
				idx=d[:,1]==t
				Rh=d[idx,6]
				n=d[idx,2]
				pl.plot(d[idx,0],1/Rh,label="T="+str(t))
				pl.plot(d[idx,0],n,label="T="+str(t))
	def getNelect(self,outcar='OUTCAR'):
		# Read number of electrons from OUTCAR.
		foundnelect = False
		try:
			for line in open(outcar, 'r'):
				if 'NELECT' in line:
					nelect = float(line.split()[2])
				foundnelect = True
		except:
			pass
		if not foundnelect:
			print 'Number of electrons not found. Please set the number manually in hte.intrans.'
			nelect = 1.0
		return nelect