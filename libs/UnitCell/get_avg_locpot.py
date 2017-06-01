#!/usr/bin/env python

from unitcell import *
import re,sys
import numpy as np

poscar = UnitCell("POSCAR")
if len(sys.argv) > 1:
    center = int(sys.argv[1])
    center = poscar.atom_positions[int(sys.argv[1])-1]
else:
    center = np.array([0.0,0.0,0.0])
if len(sys.argv) > 1+poscar.num_atom_types:
    names = [ str(i) for i in sys.argv[2:2+poscar.num_atom_types] ]
    poscar.set_atom_names(names)
if len(sys.argv) > 4+poscar.num_atom_types:
    center = [ float(i) for i in
            sys.argv[2+poscar.num_atom_types:6+poscar.num_atom_types] ]
poscar.convention = "D"
poscar.shift(0.5-center[0],0.5-center[1],0.5-center[2],"D")
poscar.in_cell()
poscar.scale = 1.0
poscar.convention = "C"



outcar = open("OUTCAR","r")
lines = outcar.readlines()
outcar.close()

avg_loc_pots = np.zeros(poscar.num_atoms)

for i in range(len(lines)):
    if lines[i].startswith("  (the norm of the test charge is"):
        j = 1
        k = 0
        while k < poscar.num_atoms:
            line = lines[i+j].split()
            while len(line) > 0:
                k = int(line.pop(0))
                avg_loc_pots[k-1] = float(line.pop(0))
            j += 1
        break
for i in range(poscar.num_atoms):
        print i+1,poscar.atom_names[i],poscar.atom_positions[i,0],poscar.atom_positions[i,1],poscar.atom_positions[i,2],avg_loc_pots[i]
sys.exit()
