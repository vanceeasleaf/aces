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
		print >>f,"velocity all create %f %d mom yes rot yes dist gaussian"%(m.T,m.seed)
		print >>f,"fix getEqu  all  nvt temp %f %f %f"%(m.T,m.T,m.dtime)
		print >>f,"run %d"%m.equTime
		print >>f,"unfix getEqu"
		print >>f,"reset_timestep 0"
		print >>f,"fix nve all nve"
		print >>f,"dump lala all custom %s velocity.txt id type vx vy vz"%m.Cinterval
		print >>f,"dump_modify  lala sort id"
		print >>f,"run %s"%m.Ctime

		f.close()
		passthru(config.mpirun+"  %s "%self.m.cores+config.lammps+" <correlation.lmp  >out.dat")
		self.vd=vdos(m.timestep)
		self.vd.run()
		#rm("velocity.txt")