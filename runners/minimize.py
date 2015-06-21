#encoding:utf8
def input(m):
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
