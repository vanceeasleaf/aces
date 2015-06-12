#encoding:utf8
import sys
from devices import nvtDevice,mpDevice,ijDevice,gkDevice
from aces.Units import Units
from aces.tools import *
import aces.config as config
from ase.io import read
from ase.io.vasp import write_vasp
def correlation(m):
	units=Units(m.units)
	m.kb=units.boltz
	m.nktv=units.nktv2p
	if(m.method=="nvt"):m.xp=0;
	lx,ly,lz,m.zfactor,m.S,xlo,xhi,ylo,yhi,zlo,zhi=postMini(m.xp,m.yp,m.zp,m.enforceThick,m.thick)
	m.dtime=m.timestep*100;
	m.tcfactor=units.tcfactor;
	m.excNum=m.aveRate/m.excRate;
	m.swapEnergyRate=m.swapEnergy/(m.excRate*m.timestep);
	f=open("correlation.lmp","w")
	print >>f,"units %s"%units
	print >>f,"dimension 3"
	pbcx=pbcy=pbcz='s'
	if m.xp==1:pbcx='p'
	if m.yp==1:pbcy='p'
	if m.zp==1:pbcz='p'
	print >>f,"boundary %s %s %s"%(pbcx,pbcy,pbcz)
	print >>f,"atom_style atomic"
	print >>f,"read_restart   minimize/restart.minimize"
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
	passthru(config.lammps+" <correlation.lmp >out.dat")
	from aces.vdos import vdos
	vdos(m.timestep)
	rm("velocity.txt")