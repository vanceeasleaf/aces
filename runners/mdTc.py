#encoding:utf8
import sys
from aces.runners.devices import nvtDevice,mpDevice,ijDevice,gkDevice
from aces.Units import Units
from aces.tools import *
import aces.config as config
from ase.io import read
from ase.io.vasp import write_vasp
from aces.runners import Runner
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
class runner(Runner):
	def generate(self):
		m=self.m
		hook=Hook()
		if m.method=='nvt':
			device=nvtDevice(hook,m)
		elif m.method=='muller':
			device=mpDevice(hook,m)
		elif m.method=='inject':
			device=ijDevice(hook,m)
		elif m.method=='greenkubo':
			device=gkDevice(hook,m)
		debug("tcfactor="+str(m.tcfactor))
		#settings
		print "units %s"%m.units
		print "dimension 3"
		pbcx=pbcy=pbcz='f'
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
		print "compute         stress all stress/atom NULL virial"
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

	def runcmd(self):
		return config.mpirun+"  %s "%self.m.cores+config.lammps+" <input  >log.out"
	


