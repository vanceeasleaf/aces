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
		cp('../minimize/POSCAR','.')
		self.m.ismear=-1
		self.getVaspRun_vasp()
		cd('..')
	def drawdos(self):
		from aces.ElectronicDOS.electronicdos import ElectronicDOS
		doscar = ElectronicDOS()
		dos=doscar.write_dos(doscar.tot_dos)
		write(dos,'dos.txt')
		dos=np.loadtxt('dos.txt');
		with fig("dos.png"):
			pl.xlabel("Energy (eV)")
			pl.ylabel("DOS")
			pl.plot(dos[:,0],dos[:,1])
	def get_fermi(self):
		a=shell_exec("grep fermi band/OUTCAR|tail -1");
		from aces.scanf import sscanf
		a= sscanf(a,"E-fermi :   %f     XC(G=0): %f     alpha+bet :%f")[0]
		print "E-fermi="+str(a)+"eV"
		return a
	def get_outfermi(self):
		file=ls("*.outputtrans")[0]
		a=shell_exec("grep FermiE %s|tail -1"%file);
		from aces.scanf import sscanf
		a= sscanf(a,"FermiE:  %f.")[0]
		print "E-fermi="+str(a)+"Ry"
		return a
	def get_nelect(self):
		a=shell_exec("grep NELECT band/OUTCAR|tail -1");
		from aces.scanf import sscanf
		a=  sscanf(a,"NELECT =      %f    total number of electrons")[0]
		print "NELECT="+str(a)
		return a
	def prepareet(self):
		fermi=self.get_fermi();
		nelect=self.get_nelect();
		mkcd('et')
		cp('../band/SYMMETRY','.')
		cp('../band/EIGENVAL','.')
		cp('../band/POSCAR','.')
		s="""VASP                      # Format of DOS                                                            
0 1 0 0.0                 # iskip (not presently used) idebug setgap shiftgap                        
%f 0.0005 0.40  %f   # Fermilevel (eV, please note), energygrid (Ry), energy span around Fermilevel (Ry), No of electrons   
CALC                      # CALC (calculate expansion coeff), NOCALC read from file                  
5                         # lpfac, number of latt-points per k-point                                      
BOLTZ                     # run mode (only BOLTZ is supported)                                       
.15                       # (efcut) energy range of chemical potential (Ry) 
300. 300.                 # Tmax, temperature grid                                                   
0.             # energyrange of bands given individual DOS output sig_xxx and dos_xxx (xxx is band number)
"""%(fermi,nelect)
		write(s,"et.intrans")
		cd('..')
	def runet(self):
		cd('et')
		m=self.m
		passthru(config.x_trans_vasp+" BoltzTrap_vasp")
		cd('..')
	def showCmd(self):
		print config.x_trans_vasp+" BoltzTrap_vasp"
	def post(self):
		file=ls("*.trace")[0]
		d=np.loadtxt(file,skiprows=1)
		head=shell_exec("head -1 %s"%file)
		if("Ry" in head):
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
		with fig("kappa.png"):#W/mK*1/s
			pl.xlabel("Energy (eV)")
			pl.ylabel("Electronic Thermal Conductivity (W/mK)")
			for t in [T[-1]]:
				idx=d[:,1]==t
				pl.plot(d[idx,0],d[idx,7]*1e-12,label="T="+str(t))
		with fig("powerfactor.png"):
			pl.xlabel("Energy (eV)")
			pl.ylabel("$S^{2}\sigma $")
			for t in [T[-1]]:
				idx=d[:,1]==t
				S=d[idx,4]
				sigma=d[idx,5]
				pl.plot(d[idx,0],S*S*sigma,label="T="+str(t))
		with fig("sigma.png"):#1/(ohm m)*1/s
			pl.xlabel("Energy (eV)")
			pl.ylabel("$\\sigma (1/\\Omega m) $")
			for t in [T[-1]]:
				idx=d[:,1]==t
				sigma=d[idx,5]*1e-12 #sigma/tau*1ps
				pl.plot(d[idx,0],sigma,label="T="+str(t))
		with fig("Rh-n.png"):
			pl.xlabel("Energy (eV)")
			pl.ylabel("$\\sigma/\\tau $")
			for t in [T[-1]]:
				idx=d[:,1]==t
				Rh=d[idx,6]
				n=d[idx,2]
				pl.plot(d[idx,0],1/Rh,label="T="+str(t))
				pl.plot(d[idx,0],n,label="T="+str(t))