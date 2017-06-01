#!/usr/bin/env python
from unitcell import *
import sys

poscar = open(sys.argv[1],'r')
unit_cell = UnitCell(poscar)
outfile = unit_cell.spacegroup()
print outfile
sys.exit()
