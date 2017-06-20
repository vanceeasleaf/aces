# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2017-06-07 14:59:50
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-20 15:15:09

import aces.tools as tl
from aces import config
import numpy as np


def minimize(m):
    if m.engine == "lammps":
        minimize_lammps(m)
    elif m.engine == "vasp":
        minimize_vasp(m)
    else:
        raise Exception('Unknow minimize type!')


def minimize_vasp(m):
    npar = 1
    for i in range(1, int(np.sqrt(m.cores)) + 1):
        if m.cores % i == 0:
            npar = i
    if m.ispin:
        ispin = "ISPIN=2"
    else:
        ispin = ""
    if m.soc:
        soc = "LSORBIT=T"
    else:
        soc = ""
    if m.isym:
        sym = "ISYM = 1"
    else:
        sym = "ISYM = 0"
    s = """SYSTEM = - local optimisation
	PREC = high
	ENCUT=%f
	EDIFF = 1e-8
	IBRION = 2
	NSW=100
	ISIF = 3
	ISMEAR = 0 ; SIGMA = 0.1
	POTIM=0.01
	ISTART = 0
	LWAVE = FALSE
	LCHARG = FALSE
	EDIFFG = -0.01
	LREAL=FALSE
	NPAR = %d
	%s
	%s
	%s
	""" % (m.ecut, npar, sym, ispin, soc)
    if m.vdw:
        s += """\nIVDW = 1
		VDW_RADIUS = 50
		VDW_S6 = 0.75
		VDW_SR = 1.00
		VDW_SCALING = 0.75
		VDW_D = 20.0
		VDW_C6 = 63.540 31.50
		VDW_R0 = 1.898 1.892
		"""
    tl.write(s, 'INCAR')
    m.structure()
    from aces.io.vasp import writePOTCAR
    writePOTCAR(m, m.elements)
    s = """A
	0
	Monkhorst-Pack
	%s
	0  0  0
	""" % ' '.join(map(str, m.mekpoints))
    tl.write(s, 'KPOINTS')
    vasp = [config.vasp, config.vasp_2d][m.d2]
    if m.useMini:
        tl.shell_exec(config.mpirun + " %s " % m.cores + vasp + ' >log.out')
    else:
        tl.cp('POSCAR', 'CONTCAR')


def minimize_lammps(m):
    f = open('input', 'w')
    print >> f, "units %s" % m.units
    print >> f, m.getatomicstyle()
    print >> f, "boundary p p p"
    print >> f, "dimension 3"
    m.structure()
    print >> f, 'read_data structure'
    print >> f, m.potential
    print >> f, "timestep %f" % m.timestep
    print >> f, m.masses
    print >> f, "thermo_style custom step pe etotal"
    print >> f, "thermo %d" % m.dumpRate
    if m.write_structure:
        print >> f, "write_data structure"
        print >> f, "dump dumpc all xyz 1CN.xyz"
        print >> f, "run 0"
    if m.metropolis:
        print >> f, "min_style metropolis"
        print >> f, "minimize 1e-12 1e-12 1000000 1000000"
    if m.useMini:
        if m.boxOpt:
            if m.enforceThick:
                print >> f, "fix 1 all box/relax x 0.0 y 0.0 nreset 1"
            else:
                print >> f, "fix 1 all box/relax iso 0.0 nreset 1"
        print >> f, "min_style cg"
        print >> f, "minimize 1e-12 1e-12 1000000 1000000"
    print >> f, "write_restart restart.minimize"
    print >> f, "dump dump1 all xyz 1 minimize.xyz"
    print >> f, m.dump
    print >> f, "dump kaka all atom 1 range"
    print >> f, "dump_modify  kaka sort id"
    print >> f, "dump 1 all custom 1 dump.force id  fx fy fz xs ys zs"
    print >> f, "dump_modify 1 format \"%d %f %f %f %f %f %f\""
    print >> f, "dump_modify  1 sort id"
    print >> f, "run 0"
    f.close()
    tl.shell_exec(config.mpirun + " %s " % m.cores + config.lammps +
                  " <input >log.out")
    m.postMini()
