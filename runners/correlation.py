#encoding:utf8

from aces.tools import *
import aces.config as config

from aces.runners.vdos import vdos
from aces.runners import Runner
from ase.io import read
from aces.lammpsdata import lammpsdata
class runner(Runner):
	def get_structure(self):
		atoms=self.m.atoms_from_dump('minimize/range')
		atoms=atoms.repeat(self.m.correlation_supercell)
		a=lammpsdata(atoms,self.m.elements)
		a.writedata('correlation_structure')
	def generate(self):
		m=self.m
		self.get_structure()
		f=open("correlation.lmp","w")
		print >>f,"units %s"%m.units
		print >>f,"dimension 3"
		pbcx=pbcy=pbcz='s'
		if m.xp==1:pbcx='p'
		if m.yp==1:pbcy='p'
		if m.zp==1:pbcz='p'
		print >>f,"boundary %s %s %s"%(pbcx,pbcy,pbcz)
		print >>f,"atom_style atomic"
		print >>f,"read_data   correlation_structure"
		print >>f,"change_box	all	boundary %s %s %s"%(pbcx,pbcy,pbcz)
		print >>f,"lattice fcc 5" #needed to define the regions
		print >>f,"thermo %d"%m.dumpRate
		print >>f,"thermo_modify     lost warn"
		print >>f,m.masses
		print >>f,m.potential
		print >>f,"timestep %f"%m.timestep
		print >>f,"reset_timestep 0"
		
		if m.corrNVT:
			box=m.box
			deta=m.deta
			wfix=m.wfix
			xlo,xhi,ylo,yhi,zlo,zhi,lx,ly,lz=box
			fixl1=xlo-deta;fixl2=fixl1+deta*wfix
			fixr2=xhi+deta;fixr1=fixr2-deta*wfix
			print >>f,"region	stayl	block   %s  %s INF  INF INF  INF units box"%(fixl1,fixl2)
			print >>f,"region	stayr	block   %s  %s INF INF   INF  INF units box"%(fixr1,fixr2)
			print >>f,"region   stay    union  2 stayl stayr"
			print >>f,"region	main	block   %s  %s INF INF   INF  INF units box"%(fixl2,fixr1)
			print >>f,"group   stay    region  stay"
			print >>f,"group   main    region  main"
			print >>f,"velocity stay set 0 0 0"
		else:
			print >>f,"region	main	block   INF  INF INF  INF INF  INF units box"
			print >>f,"group    main    region  main"
		print >>f,"velocity main create %f %d mom yes rot yes dist gaussian"%(m.T,m.seed)
		print >>f,"fix getEqu  main  nvt temp %f %f %f"%(m.T,m.T,m.dtime)
		print >>f,"dump dump1 all atom %d dump.lammpstrj"%(m.Cinterval*500)
		print >>f,"dump_modify  dump1 sort id"
		print >>f,"run %d"%m.equTime
		print >>f,"unfix getEqu"
		print >>f,"reset_timestep 0"
		print >>f,"fix nve main nve"
		print >>f,"dump lala main custom %s velocity.txt id type vx vy vz"%m.Cinterval
		print >>f,"dump_modify  lala sort id"

		print >>f,"run %s"%m.Ctime

		f.close()
		passthru(config.mpirun+"  %s "%self.m.cores+config.lammps+" <correlation.lmp  >out.dat")
		self.vd=self.getvdos()
		self.vd.run()
		#rm("velocity.txt")
	def getvdos(self):
		return vdos(self.m.timestep)