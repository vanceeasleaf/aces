#encoding:utf8
import sys
from devices import nvtDevice,mpDevice,ijDevice,gkDevice
from aces.Units import Units
from aces.tools import *
import aces.config as config
from ase.io import read
from ase.io.vasp import write_vasp
def getboxrange():
	file=open("minimize/range");
	for i in range(5):
		file.next()
	xlo,xhi=map(float,file.next().split()[:2])
	ylo,yhi=map(float,file.next().split()[:2])
	zlo,zhi=map(float,file.next().split()[:2])
	return (xlo,xhi,ylo,yhi,zlo,zhi);

def getxrange():
    file=open('minimize/minimize.xyz');
    n=file.next().split()[0];n=int(n)
    file.next()
    xmin=100000;xmax=-100000;
    ymin=100000;ymax=-100000;
    zmin=100000;zmax=-100000;
    for i in range(n):
        label,x,y,z=file.next().split();
        x,y,z=map(float,[x,y,z]);
        xmin=min(x,xmin);xmax=max(x,xmax);
        ymin=min(y,ymin);ymax=max(y,ymax);
        zmin=min(z,zmin);zmax=max(z,zmax); 
    return (xmin,xmax,ymin,ymax,zmin,zmax);
    
def postMini(xp,yp,zp,enforceThick,thick):
	xlo,xhi,ylo,yhi,zlo,zhi=getboxrange();
	xlo0,xhi0,ylo0,yhi0,zlo0,zhi0=getxrange();
	if(xp==0):
		xlo=xlo0;xhi=xhi0;
	if(yp==0):
		ylo=ylo0;yhi=yhi0;
	if(zp==0):
		zlo=zlo0;zhi=zhi0;
	lx=xhi-xlo;ly=yhi-ylo;lz=zhi-zlo;
	if(enforceThick):zfactor=lz/thick;
	else:zfactor=1;
	S=ly*lz;
	return lx,ly,lz,zfactor,S,xlo,xhi,ylo,yhi,zlo,zhi



class Hook:
	def __init__(self):
		self.labels={}
	def addAction(self,label,function):
		if not self.labels.has_key(label):
			self.labels[label]=[]
		self.labels[label].append(function)
	def doAction(self,label):
		if not self.labels.has_key(label):return
		for key in self.labels[label]:
			key()
def mdTc(m):
	units=Units(m.units)
	m.kb=units.boltz
	m.nktv=units.nktv2p
	if(m.method=="nvt"):m.xp=0;
	lx,ly,lz,m.zfactor,m.S,xlo,xhi,ylo,yhi,zlo,zhi=postMini(m.xp,m.yp,m.zp,m.enforceThick,m.thick)
	m.dtime=m.timestep*100;
	m.tcfactor=units.tcfactor;
	m.excNum=m.aveRate/m.excRate;
	m.swapEnergyRate=m.swapEnergy/(m.excRate*m.timestep);

	m.box=(xlo,xhi,ylo,yhi,zlo,zhi,lx,ly,lz)
	hook=Hook()
	if m.method=='nvt':
		device=nvtDevice(hook,m)
	elif m.method=='muller':
		device=mpDevice(hook,m)
	elif m.method=='inject':
		device=ijDevice(hook,m)
	elif m.method=='greenkubo':
		m.corNum=m.aveRate/m.corRate;
		device=gkDevice(hook,m)

	#settings
	print "units %s"%units
	print "dimension 3"
	pbcx=pbcy=pbcz='s'
	if m.xp==1:pbcx='p'
	if m.yp==1:pbcy='p'
	if m.zp==1:pbcz='p'
	print "boundary %s %s %s"%(pbcx,pbcy,pbcz)
	print "atom_style atomic"
	print "read_restart   minimize/restart.minimize"
	print "change_box	all	boundary %s %s %s"%(pbcx,pbcy,pbcz)
	print "lattice fcc 5" #needed to define the regions
	print "thermo %d"%m.dumpRate
	print "thermo_modify     lost warn"
	print "timestep %f"%m.timestep
	#regions and groups
	hook.doAction('region')
	#computes
	print "compute           ke  all  ke/atom"
	print "compute           pe  all  pe/atom"
	print "compute         stress all stress/atom virial"
	print "compute jflux all heat/flux ke pe stress"
	hook.doAction('compute')
	#variables
	hook.doAction('variable')
	#init atoms to T
	print m.masses
	print m.potential
	print "reset_timestep 0"
	print "velocity all create %f %d mom yes rot yes dist gaussian"%(m.T,m.seed)
	hook.doAction('equ')
	print "run %d"%m.equTime
	print "unfix getEqu"
	print "reset_timestep 0"
	hook.doAction('elimination')
	print "fix    flux_out  all  ave/time  1  %d  %d  c_jflux[1]  c_jflux[2] c_jflux[3] file  flux.txt "%(m.aveRate,m.aveRate)
	hook.doAction('temp')
	hook.doAction('flux')
	hook.doAction('swap')


	#/* 定时输出dump文件并按id排序*/
	if(m.dumpxyz):
		print "dump dump1 all atom %d dump.lammpstrj"%(m.dumpRate)
		print "dump_modify  dump1 sort id"


	#/* 定时输出速度文件用于计算速度关联函数*/
	if(m.dumpv):
		print "dump dump2 all custom %d dump.velocity type vx vy vz"%(m.dumpRate)
		print "dump_modify  dump2 sort id"

	print "run	%d"%(m.runTime)

def input(m):
	if m.bte==True:
		from aces.bte import bte
		bte(m)
	elif m.bte=="phonopy":
		from aces.phonopy import phonopy
		phonopy(m)
	elif m.bte=="shengbte":shengbte(m)
	elif m.bte=="correlation":
		from aces.correlation import correlation
		correlation(m)
	else:mdTc(m)

def shengbte(m):
	pass	


