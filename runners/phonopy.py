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
class runner(Runner):
	def minimizePOSCAR(self):
		m=self.m
		if m.engine=="lammps":

			m.dump2POSCAR('minimize/range')
			
		elif m.engine=="vasp":
			cp('minimize/CONTCAR','POSCAR')
	
	def get_lammps_script(self,m):
		content="units %s\n"%m.units
		content+="""atom_style      atomic
dimension       3
boundary        p p p 
read_data       structure
%s
%s 
neighbor        2 bin
neigh_modify    every 1 delay 1 check yes
dump 1 all custom 1 dump.force id  fx fy fz xs ys zs
dump_modify 1 format "%%d %%30.20f %%30.20f %%30.20f %%30.20f %%30.20f %%30.20f"
dump_modify  1 sort id
run 0
"""%(m.masses,m.potential)
		return content
		
	def force_constant(self,files):
		cmd=config.phonopy+"-f "
		for file in files:
			dir="dirs/dir_"+file
			cmd+=dir+'/vasprun.xml '
		#generate FORCE_SETS
		passthru(cmd)
		
	def generate_meshconf(self):
		#generate mesh.conf
		m=self.m

		mesh="""DIM = %s
ATOM_NAME = %s
MP = %s
EIGENVECTORS=.TRUE.
FORCE_CONSTANTS = WRITE
"""%(m.dim,' '.join(m.elements),' '.join(map(str,m.kpoints)))
		write(mesh,'mesh.conf')
	
	def generate_supercells(self):
		m=self.m
		#generate supercells

		passthru(config.phonopy+"-d --dim='%s'"%(m.dim))
	
	def getVaspRun_lammps(self):
		m=self.m
		
		#generate structure
		m.POSCAR2data()
		#generate in
		content=self.get_lammps_script(m)
		write(content,"in")
		#generate dump.force
		shell_exec(config.lammps+" < in >log.out")
		#generate vasprun.xml
		f=open('dump.force')
		for i in range(9):f.next()
		forces=""
		poses=""
		for line in f:
			line=line.split()
			forces+="<v>  %s %s %s </v>\n"%tuple(line[1:4])
			poses+="<v>  %s %s %s  </v>\n"%tuple(line[4:8])
		vasprun='<root><calculation><varray name="forces" >\n'
		vasprun+=forces
		vasprun+='</varray>\n<structure><varray name="positions">\n'+poses
		vasprun+='</varray></structure></calculation></root>\n'
		write(vasprun,'vasprun.xml')
		
	def getVaspRun_vasp(self):
		s="""SYSTEM=calculate energy
PREC = Accurate
IBRION = -1
ENCUT = 300
EDIFF = 1.0e-6
ISMEAR = 0; SIGMA = 0.01
IALGO = 38
LREAL = .FALSE.
ADDGRID = .TRUE.
LWAVE = .FALSE.
LCHARG = .FALSE.
"""
		write(s,'INCAR')
		m=self.m
		m.writePOTCAR()
		s="""A
0
Monkhorst-Pack
%s
0  0  0
	"""%' '.join(map(str,m.kpoints))
		write(s,'KPOINTS')
		shell_exec(config.mpirun+" %s "%m.cores+config.vasp+' >log.out')
	
	def getvasprun(self,files):
		m=self.m
		maindir=shell_exec('pwd')
		for file in files:
			print file
			dir="dirs/dir_"+file
			mkdir(dir)
			mv(file,dir+'/POSCAR')
			cd(dir)
			if m.engine=="lammps":			
				self.getVaspRun_lammps()
			elif m.engine=="vasp":
				self.getVaspRun_vasp()
			cd(maindir)

	def generate(self):
		m=self.m
		self.minimizePOSCAR()
		self.generate_supercells()
		files=shell_exec("ls *-*").split('\n')
		self.getvasprun(files)
		self.force_constant(files)
		
		self.generate_meshconf()	
		passthru(config.phonopy+" --dos  mesh.conf")
		self.drawDos()
		
		self.generate_bandconf()
		passthru(config.phonopy+" -s  band.conf")
		from aces.bandplot import plotband
		plotband(labels=' '.join(m.bandpath))
	def generate_bandconf(self):
		#generate mesh.conf
		m=self.m

		bp=m.bandpoints
		bpath=' '.join([' '.join(map(str,bp[x])) for x in m.bandpath])
		
		band="""DIM = %s
ATOM_NAME = %s
BAND = %s 
BAND_POINTS = 101
"""%(m.dim,' '.join(m.elements),bpath)
		write(band,'band.conf')
		
	def drawDos(self):
		xx=np.loadtxt('partial_dos.dat',skiprows=1)
		ndos=len(xx[0,:])-1
		freq=xx[:,0]
		pdos=xx[:,1:]
		ndos=len(pdos[0,:])
		
		datas=[]
		for i in range(ndos):
			datas.append((freq,pdos[:,i],''))
		series('Frequency (THz)','Partial Density of States',
		datas=datas,
		filename='partial_dos.png',legend=False,grid=True)
		
		plot((freq,'Frequency (THz)'),(np.sum(pdos,axis=1),'Density of States'),filename='total_dos.png')
		#calculate paticipation ratio
		
		pr()
		#plot
		xs=[];ys=[]
		for line in open('pr.txt'):
			x,y=map(float,line.split())
			xs.append(x);ys.append(y)
		write("%s"%(sum(ys)/len(ys)),"ave_pr.txt")
		pl.figure()
		pl.plot(xs,ys,'.',color='red')
		pl.ylim([0.0,1.0])
		pl.xlabel('Frequency (THz)')
		pl.ylabel('Paticipation Ratio')
		pl.savefig('Paticipation_atio.png',bbox_inches='tight',transparent=True) 
		pl.close()