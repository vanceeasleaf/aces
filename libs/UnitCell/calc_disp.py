#!/usr/bin/env python

from unitcell import *
import sys

cell_1 = UnitCell(open(str(sys.argv[1]),'r'))
cell_2 = UnitCell(open(str(sys.argv[2]),'r'))

#disp = displacements(cell_1,cell_2,conv='C')
#disp = np.round(disp,decimals=9)
mag,dir = displacements(cell_1,cell_2,conv='C',flag=True)
print mag
