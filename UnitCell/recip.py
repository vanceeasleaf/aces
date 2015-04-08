#!/usr/bin/python
from unitcell import *
import sys

file = sys.argv[1]
cell_file = open(file,'r')
cell = UnitCell(cell_file)
recip,b = cell.recip_lat()
print recip
print
print b
