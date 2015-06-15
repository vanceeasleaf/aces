#encoding:utf8
def input(m):
	units,structure,potential,timestep,masses,dumpRate,write_structure,metropolis,useMini,dump=m.units,m.structure,m.potential,m.timestep,m.masses,m.dumpRate,m.write_structure,m.metropolis,m.useMini,m.dump
	print "units %s"%units
	print "atom_style atomic"
	print "boundary p p p"
	print "dimension 3"
	structure()
	print potential
	print "timestep %f"%timestep
	print masses
	print "thermo_style custom step pe etotal"
	print "thermo %d"%dumpRate
	if write_structure:
		print "write_data structure"
		print "dump dumpc all xyz 1CN.xyz"
		print "run 0"
	if metropolis:
		print "min_style metropolis"
		print "minimize 1e-12 1e-12 1000000 1000000"
	if useMini:
		print "fix 1 all box/relax x 0.0 y 0.0 nreset 1"
		print "min_style cg"
		print "minimize 1e-12 1e-12 1000000 1000000"
	print "write_restart restart.minimize"
	print "dump dump1 all xyz 1 minimize.xyz"
	print dump
	print "dump kaka all atom 1 range"
	print "run 0"
