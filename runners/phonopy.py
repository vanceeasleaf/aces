#encoding:utf8
from aces.tools import *
import aces.config as config
from ase.io import read
from ase.io.vasp import write_vasp
from aces.binary import pr
from aces.runners import Runner
from aces.UnitCell.unitcell import UnitCell
from aces.graph import plot,series
from aces.script.vasprun import exe as lammpsvasprun
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pl
import time
import numpy as np

class runner(Runner):
	def minimizePOSCAR(self):
		m=self.m
		if m.engine=="lammps":

			m.dump2POSCAR(m.home+'/minimize/range',rotate=True)
			
		elif m.engine=="vasp":
			cp(m.home+'/minimize/CONTCAR','POSCAR')
	
		
	def force_constant(self,files):
		cmd=config.phonopy+"-f "
		for file in files:
			dir="dirs/dir_"+file
			cmd+=dir+'/vasprun.xml '
		#generate FORCE_SETS
		passthru(cmd)
		m=self.m
		#Create FORCE_CONSTANTS
		passthru(config.phonopy+"--tolerance=1e-4 --writefc --dim='%s'"%(m.dim))

	def generate_meshconf(self):
		#generate mesh.conf
		m=self.m

		mesh="""DIM = %s
ATOM_NAME = %s
MP = %s
EIGENVECTORS=.TRUE.
FORCE_CONSTANTS = WRITE
MESH_SYMMETRY = .FALSE.
PRIMITIVE_AXIS = %s
"""%(m.dim,' '.join(m.elements),' '.join(map(str,m.kpoints)),toString(m.premitive.flatten()))
		write(mesh,'mesh.conf')
	def generate_vconf(self):
		#generate v.conf
		m=self.m

		mesh="""DIM = %s
ATOM_NAME = %s
MP = %s
FORCE_CONSTANTS = READ
MESH_SYMMETRY = .FALSE.
GROUP_VELOCITY=.TRUE.
PRIMITIVE_AXIS = %s
"""%(m.dim,' '.join(m.elements),' '.join(map(str,m.kpoints)),toString(m.premitive.flatten()))
		write(mesh,'v.conf')	
	def generate_qconf(self,q):
		#generate q.conf
		m=self.m

		mesh="""DIM = %s
ATOM_NAME = %s
FORCE_CONSTANTS = READ
EIGENVECTORS=.TRUE.
QPOINTS=.TRUE.
PRIMITIVE_AXIS = %s
"""%(m.dim,' '.join(m.elements),toString(m.premitive.flatten()))
		write(mesh,'q.conf')	
		s="%s\n"%len(q)
		for qq in q:
			s+="%s\n"%m.toString(qq)
		write(s,'QPOINTS')
	def generate_vqconf(self,q):
		#generate q.conf
		m=self.m

		mesh="""DIM = %s
ATOM_NAME = %s
FORCE_CONSTANTS = READ
GROUP_VELOCITY=.TRUE.
QPOINTS=.TRUE.
PRIMITIVE_AXIS = %s
"""%(m.dim,' '.join(m.elements),toString(m.premitive.flatten()))
		write(mesh,'q.conf')	
		s="%s\n"%len(q)
		for qq in q:
			s+="%s\n"%m.toString(qq)
		write(s,'QPOINTS')
	def generate_supercells(self):
		m=self.m
		#generate supercells

		passthru(config.phonopy+"--tolerance=1e-4 -d --dim='%s'"%(m.dim))

	
		
	def getVaspRun_vasp(self):
		m=self.m 
		if m.isym:
			sym="ISYM = 1"
		else:
			sym="ISYM = 0"
		s="""SYSTEM=calculate energy
PREC = High
IBRION = -1
ENCUT = %f
EDIFF = 1.0e-8
ISMEAR = -5; SIGMA = 0.01
IALGO = 38
LREAL = .FALSE.
ADDGRID = .TRUE.
LWAVE = .FALSE.
LCHARG = .FALSE.
NPAR = %d
%s
"""%(self.m.ecut,m.nodes,sym)
		write(s,'INCAR')
		m=self.m
		m.writePOTCAR()
		s="""A
0
Monkhorst-Pack
%s
0  0  0
	"""%' '.join(map(str,m.ekpoints))
		write(s,'KPOINTS')
		if 'jm' in self.__dict__:
			from aces.jobManager import pbs
			path=pwd()
			if m.queue=="q3.4":
				pb=pbs(queue=m.queue,nodes=12,procs=1,disp=m.pbsname,path=path,content=config.mpirun+" 12 "+config.vasp+' >log.out')
			else:
				pb=pbs(queue=m.queue,nodes=1,procs=12,disp=m.pbsname,path=path,content=config.mpirun+" 12 "+config.vasp+' >log.out')
			self.jm.reg(pb)
		else:
			shell_exec(config.mpirun+" %s "%m.cores+config.vasp+' >log.out')
	
	def getvasprun(self,files):
		m=self.m
		maindir=pwd()
		if m.engine=="vasp":
			from aces.jobManager import jobManager
			self.jm=jobManager()
			for file in files:
				print file
				dir="dirs/dir_"+file
				mkdir(dir)
				mv(file,dir+'/POSCAR')
				cd(dir)
				self.getVaspRun_vasp()
				cd(maindir)
			self.jm.check()
		elif m.engine=="lammps":
			from multiprocessing.dummy  import Pool
			pool=Pool()
			pool.map_async(lammpsvasprun,files)
			pool.close()
			pool.join()



	def runSPOSCAR(self):
		m=self.m
		maindir=pwd()
		file="SPOSCAR"
		dir="dir_"+file
		mkdir(dir)
		cp(file,dir+'/POSCAR')
		cd(dir)		
		self.getVaspRun_lammps()
		cd(maindir)
	def checkMinimize(self):
		import yaml
		data=yaml.load(open('disp.yaml').read())
		disps=[map(float,a['direction']) for a in data['displacements']]
		maindir=pwd()
		dirs=ls('dirs/dir_*')
		ii=0
		L=np.linalg.norm
		#d,p,d1,p1=self.m.rot
		out=open('ccos.txt','w')
		for dir in dirs:
			cd(dir)
			f=open('dump.force')
			for i in range(9):f.next()
			forces=""
			for b in range(ii):f.next()
			line=f.next()
			line=line.split()
			force=np.array(map(float,line[1:4]))
			#force=RotateVector(force,d1,-p1)
			#force=RotateVector(force,d,-p)
			d=disps[i]
			ccos=force.dot(d)/L(force)/L(d)
			ii+=1
			print >>out,"%d\t%f"%(ii,ccos)
			cd(maindir)
	def generate(self):
		m=self.m
		self.minimizePOSCAR()
		a=time.time()
		self.generate_supercells()
		debug('generate_supercells:%f s'%(time.time()-a))
		files=shell_exec("ls *-*").split('\n')
		assert len(files)>0
		#self.runSPOSCAR()
		a=time.time()
		self.getvasprun(files)
		debug('getvasprun:%f s'%(time.time()-a))
		a=time.time()
		self.force_constant(files)
		debug('force_constant:%f s'%(time.time()-a))
		
		if m.phofc:return self
		self.postp()
	def get_force_sets(self):
		files=shell_exec("ls dirs").split('\n')
		files=map(lambda x:x.replace('dir_',''),files)
		self.force_constant(files)
	def postp(self):
		self.getband()
		self.getDos()
		
		
		self.getbanddos()
		self.drawpr()
		self.getV()
	def getqpoints(self,q):
		self.generate_qconf(q)
		passthru(config.phonopy+"--tolerance=1e-4 q.conf")
	def getvqpoints(self,q):
		self.generate_vqconf(q)
		passthru(config.phonopy+"--tolerance=1e-4 q.conf")
		data =parseyaml('qpoints.yaml')
		file=open("v.txt",'w')
		for phonon in data['phonon']:
			qp=phonon['q-position']
			for band in phonon['band']:
				frequency=band['frequency']
				v=np.array(band['group_velocity'])
				v=np.linalg.norm(v)
				print >>file,"%s\t%f\t%f"%('\t'.join(map(str,qp)),frequency,v)
		file.close()
		v=np.loadtxt('v.txt')
		plot((v[:,3],'Frequency (THz)'),(v[:,4],'Group Velocity (Angstrom/ps)'),'v_freq.png',grid=True,scatter=True)
	def getDos(self):
		self.generate_meshconf()	
		passthru(config.phonopy+"--tolerance=1e-4 --dos  mesh.conf")
		self.drawDos()
	def getV(self):
		if not exists('groupv'):mkdir('groupv')
		cd('groupv')
		cp('../FORCE_CONSTANTS','.')
		cp('../POSCAR','.')
		cp('../disp.yaml','.')
		self.generate_vconf()	
		passthru(config.phonopy+"--tolerance=1e-4   v.conf")
		self.drawV()
		cd('..')
	def drawV(self):
		data =parseyaml('mesh.yaml')
		file=open("v.txt",'w')
		for phonon in data['phonon']:
			qp=phonon['q-position']
			for band in phonon['band']:
				frequency=band['frequency']
				v=np.array(band['group_velocity'])
				v=np.linalg.norm(v)
				print >>file,"%s\t%f\t%f"%('\t'.join(map(str,qp)),frequency,v)
		file.close()
		v=np.loadtxt('v.txt')
		plot((v[:,3],'Frequency (THz)'),(v[:,4],'Group Velocity (Angstrom/ps)'),'v_freq.png',grid=True,scatter=True)
	def getbandpr(self):
		pr=np.loadtxt('pr.txt')
		v=np.loadtxt('groupv/v.txt')
		
	def getband(self):
		self.generate_bandconf()
		passthru(config.phonopy+"--tolerance=1e-4 -s  band.conf")
		from aces.bandplot import plotband
		plotband(labels=' '.join(self.m.bandpath))
	def getbanddos(self):
		freq,pdos=self.getpdos()
		from aces.bandplot import plotbanddos
		plotbanddos(freq=freq,dos=np.sum(pdos,axis=1),labels=' '.join(self.m.bandpath))
	def generate_bandconf(self):
		#generate mesh.conf
		m=self.m

		bp=m.bandpoints
		bpath=' '.join([m.toString(bp[x]) for x in m.bandpath])
		
		band="""DIM = %s
ATOM_NAME = %s
BAND = %s 
BAND_POINTS = 101
MESH_SYMMETRY = .FALSE.
PRIMITIVE_AXIS = %s
"""%(m.dim,' '.join(m.elements),bpath,toString(m.premitive.flatten()))
		write(band,'band.conf')
	def getpdos(self):
		xx=np.loadtxt('partial_dos.dat',skiprows=1)
		freq=xx[:,0]
		pdos=xx[:,1:]
		return freq,pdos
	def drawDos(self):
		freq,pdos=self.getpdos()
		datas=[(freq,p,'') for p in pdos.T]
		series('Frequency (THz)','Partial Density of States',
		datas=datas,
		filename='partial_dos.png',legend=False,grid=True)
		
		plot((freq,'Frequency (THz)'),(np.sum(pdos,axis=1),'Density of States'),filename='total_dos.png')
		#calculate paticipation ratio
	def drawpr(self):		
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
