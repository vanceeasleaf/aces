#!/usr/bin/env python
# v1.0 Jeff Doak jwd686@u.northwestern.edu 8-1-2011
# Note: Currently only supports orthogonal cells with vec_a1 being along x,
# vec_a2 being along y, and vec_a3 being along z.
from unitcell import *
import sys
poscar = open(sys.argv[1],'r')
unit_cell = UnitCell(poscar)
text = unit_cell.output_lammps()
print text
