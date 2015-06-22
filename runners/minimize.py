#encoding:utf8
from aces.tools import *
from aces import config
def minimize(m):
	if m.engine=="lammps":minimize_lammps(m)
	elif m.engine=="vasp":minimize_vasp(m)
	else: raise Exception('Unknow minimize type!')
def minimize_vasp(m):
	s="""SYSTEM = - local optimisation
PREC = accurate
ENCUT=400
EDIFF = 1e-6
IBRION = -1
ISIF = 2
ISMEAR = 0 ; SIGMA = 0.1
ISTART = 0
LWAVE = FALSE
LCHARG = FALSE
EDIFFG = 0.01
ALGO=fast
LREAL=AUTO
LPLANE=.TRUE.
"""
	write(s,'INCAR')
	m.structure()
	m.writePOTCAR()
	s="""A
0
Monkhorst-Pack
%s
0  0  0
	"""%' '.join(map(str,m.kpoints))
	write(s,'KPOINTS')
	shell_exec(config.mpirun+" %s "%m.cores+config.vasp+' >log.out')
def minimize_lammps(m):
	f=open('input', 'w')
	units,structure,potential,timestep,masses,dumpRate,write_structure,metropolis,useMini,dump=m.units,m.structure,m.potential,m.timestep,m.masses,m.dumpRate,m.write_structure,m.metropolis,m.useMini,m.dump
	print >>f,"units %s"%units
	print >>f,"atom_style atomic"
	print >>f,"boundary p p p"
	print >>f,"dimension 3"
	structure()
	print >>f,'read_data structure'
	print >>f,potential
	print >>f,"timestep %f"%timestep
	print >>f,masses
	print >>f,"thermo_style custom step pe etotal"
	print >>f,"thermo %d"%dumpRate
	if write_structure:
		print >>f,"write_data structure"
		print >>f,"dump dumpc all xyz 1CN.xyz"
		print >>f,"run 0"
	if metropolis:
		print >>f,"min_style metropolis"
		print >>f,"minimize 1e-12 1e-12 1000000 1000000"
	if useMini:
		print >>f,"fix 1 all box/relax x 0.0 y 0.0 nreset 1"
		print >>f,"min_style cg"
		print >>f,"minimize 1e-12 1e-12 1000000 1000000"
	print >>f,"write_restart restart.minimize"
	print >>f,"dump dump1 all xyz 1 minimize.xyz"
	print >>f,dump
	print >>f,"dump kaka all atom 1 range"
	print >>f,"run 0"
	f.close()
	shell_exec(config.mpirun+" %s "%m.cores+config.lammps+" <input >log.out")	
	m.postMini()
