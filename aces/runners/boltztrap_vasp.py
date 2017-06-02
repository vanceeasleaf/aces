#encoding:utf8
from aces.runners import Runner
from ase import io
from aces.tools import *
import numpy as np
from  aces import config
from aces.runners.phonopy import runner as Runner
from aces.f import read_forces,matrixFormat
import  libs.vasp2boltz
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
		#orbital_dos = doscar.sum_ms_dos()
		## Create a list of each atom type to sum over.
		#type_list = []
		#n = 0
		#for i in range(len(doscar.unit_cell.atom_types)):
		#    type_list.append([])
		#    for j in range(doscar.unit_cell.atom_types[i]):
		#        type_list[i].append(n)
		#        n += 1
		## Sum dos over sets of atoms.
		#partial_dos = doscar.sum_site_dos(type_list,orbital_dos)
		dos=doscar.write_dos([doscar.tot_dos])
		write(dos,'dos.txt')
		dos=np.loadtxt('dos.txt');
		f=shell_exec("grep fermi OUTCAR|tail -1")
		from aces.scanf import sscanf
		f=sscanf(f,"E-fermi :   %f     XC(G=0):");
		with fig("dos.png"):
			pl.xlabel("Energy-Ef (eV)")
			pl.ylabel("DOS")
			pl.plot(dos[:,0]-f,dos[:,1],lw=2)
			pl.xlim([-4,4])
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
		lpfac=5
		if(self.m.soc):lpfac=20
		s="""VASP                      # Format of DOS                                                            
0 1 0 0.0                 # iskip (not presently used) idebug setgap shiftgap                        
%f 0.0005 0.40  %f   # Fermilevel (eV, please note), energygrid (Ry), energy span around Fermilevel (Ry), No of electrons   
CALC                      # CALC (calculate expansion coeff), NOCALC read from file                  
%d                         # lpfac, number of latt-points per k-point                                      
BOLTZ                     # run mode (only BOLTZ is supported)                                       
.15                       # (efcut) energy range of chemical potential (Ry) 
800. 50.                 # Tmax, temperature grid                                                   
0.             # energyrange of bands given individual DOS output sig_xxx and dos_xxx (xxx is band number)
"""%(fermi,nelect,lpfac)
		write(s,"et.intrans")
		cd('..')
	def runet(self):
		cd('et')
		m=self.m
		if(m.soc):
			passthru(config.x_trans_vasp+" BoltzTrap_vasp  -so")
		else:
			passthru(config.x_trans_vasp+" BoltzTrap_vasp")
		cd('..')
	def showCmd(self):
		print config.x_trans_vasp+" BoltzTrap_vasp"
	def kappat(self):
		a=np.loadtxt("BTE.KappaTensorVsT_CONV")
		import matplotlib as mpl 
		mpl.rcParams['axes.color_cycle']=['#e24a33','#2A749A','#988ed5']
		with fig('T_kappa.png',legend=True):
			ts=a[:,0]
			fil=ts<=800
			ts=a[fil,0]
			k1=1.0/3*(a[fil,1]+a[fil,5]+a[fil,9])
			pl.plot(ts,k1,lw=3,label="Iso-Phonon")
			pl.xlabel("Tempeature (K)")
			pl.ylabel('Thermal Conductivity (W/mK)')
			file=ls("*.trace")[0]
			d=np.loadtxt(file,skiprows=1)
			idx=d[:,0]==d[np.abs(np.unique(d[:,0])).argmin(),0]
			tao=2.93e-14
			pl.plot(d[idx,1],d[idx,7]*tao,lw=3,label="Iso-Electron")
			pl.xlim([200,800])
	def post(self):
		import matplotlib as mpl 
		mpl.rcParams['axes.color_cycle']=['#e24a33','#2A749A','#988ed5']
		
		file=ls("*.trace")[0]
		d=np.loadtxt(file,skiprows=1)
		head=shell_exec("head -1 %s"%file)
		if("Ry" in head):
			d[:,0]-=self.get_outfermi()
			d[:,0]*=13.6
		T=np.unique(d[:,1])
		zz=d[np.abs(np.unique(d[:,0])).argmin(),0]
		Tplot=[200,300,700]
		tao=2.93e-14
		#tao=2.93e-13
		#T=[300]
		with fig("Seebeck.png",legend=True):
			#pl.style.use('ggplot')
			pl.xlabel("$\\mu$ (eV)")
			pl.ylabel("Seebeck Coefficient ($\mu$V/K)")
			pl.xlim([-1.5,1.5])
			for t in Tplot:
				idx=d[:,1]==t
				pl.plot(d[idx,0],d[idx,4],lw=3,label="T="+str(t)+"K")
			
		with fig("kappa.png",legend=True):#W/mK*1/s
			pl.xlabel("$\\mu$ (eV)")
			pl.ylabel("Electronic Thermal Conductivity (W/mK)")
			for t in Tplot:
				idx=(d[:,1]==t) * ( d[:,0]<=1.5) * ( d[:,0]>=-1.5)
				pl.plot(d[idx,0],d[idx,7]*tao,lw=3,label="T="+str(t)+"K")
			pl.xlim([-1.5,1.5])
		with fig("kappa_t.png"):
			pl.xlabel("Temperature (K)")
			pl.ylabel("Electronic Thermal Conductivity (W/mK)")
			idx=d[:,0]==zz
			pl.plot(d[idx,1],d[idx,7]*tao,lw=3)
		with fig("powerfactor.png",legend=True):
			pl.xlabel("$\\mu$ (eV)")
			pl.ylabel("$S^{2}\sigma (mW/mK^2)$")
			pl.xlim([-1.5,1.5])
			for t in Tplot:
				idx=d[:,1]==t
				S=d[idx,4]*1e-6
				sigma=d[idx,5]*tao
				pl.plot(d[idx,0],S*S*sigma*1e3,lw=3,label="T="+str(t)+"K")
		try:
			a=np.loadtxt("BTE.KappaTensorVsT_CONV")	
			k1=1.0/3*(a[:,1]+a[:,5]+a[:,9])	
			with fig("ZT.png",legend=True,ncol=1):
				pl.xlabel("$\\mu$ (eV)")
				pl.ylabel("$ZT$")
				pl.xlim([-1.5,1.5])
				pl.ylim([0,1.0])
				for t in Tplot:
					idx=d[:,1]==t
					fil=a[:,0]==t
					tc=k1[fil][0]
					S=d[idx,4]*1e-6
					sigma=d[idx,5]*tao
					po=S*S*sigma
					ke=d[idx,7]*tao
					pl.plot(d[idx,0],po*t/(ke+tc),lw=3,label="T="+str(t)+"K")
		except Exception as e:
			print e
		with fig("sigma.png",legend=True):#1/(ohm m)*1/s
			pl.xlabel("$\\mu$ (eV)")
			pl.ylabel("$\\sigma (10^6/\\Omega m) $")
			pl.xlim([-1.5,1.5])
			for t in Tplot:
				idx=d[:,1]==t
				sigma=d[idx,5]*tao 
				pl.plot(d[idx,0],sigma*1e-6,lw=3,label="T="+str(t)+"K")
				
		with fig("Rh-n.png",legend=True):
			pl.xlabel("$\\mu$ (eV)")
			pl.ylabel("$\\sigma/\\tau $")
			pl.xlim([-1.5,1.5])
			for t in Tplot:
				idx=d[:,1]==t
				Rh=d[idx,6]
				n=d[idx,2]
				pl.plot(d[idx,0],1/Rh,lw=3,label="Rh,T="+str(t)+"K")
				pl.plot(d[idx,0],n,lw=3,label="nT="+str(t)+"K")