#!/usr/bin/env python
from unitcell import *
import sys

file = sys.argv[1]
cell_file = open(file,'r')
cell = UnitCell(cell_file)
findsym = cell.output_findsym()
print findsym
