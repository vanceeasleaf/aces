#!/usr/bin/env python
from unitcell import *
import sys

cell = UnitCell(open(str(sys.argv[1]),"r"))
cell.convention = "C"
#cell.vel_convention = "C"
init_pos = cell.atom_positions
vel = cell.atom_velocities
NSW = 5.0
POTIM = 1.0

delta_r = vel*NSW*POTIM
new_pos = init_pos+delta_r
print init_pos
print
print new_pos
print
print cell.vel_convention
