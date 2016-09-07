#!/usr/bin/env python

# bands.py v1.0 9-18-2012 Jeff Doak jeff.w.doak@gmail.com

from electronicdos import *
import sys

# Read in arguments
dosdir = str(sys.argv[1])
if len(sys.argv) > 2:
    shift = float(sys.argv[2])
else:
    shift = 0.0
if len(sys.argv) > 3:
    offset = int(sys.argv[3])
else:
    offset = 0
if len(sys.argv) > 4:
    tol = float(sys.argv[4])
else:
    tol = 1e-5
# Get Bandstructure
dos = ElectronicDOS(dosdir+"/DOSCAR",dosdir+"/OUTCAR",dosdir+"/POSCAR")
dos.shift_energy(shift)
gaps = dos.get_band_gaps(tol=tol)
for i in range(len(gaps)):
    for j in range(len(gaps[i])):
        for k in range(len(gaps[i][j])):
            print i+offset,gaps[i][j][k]
            print i+1+offset,gaps[i][j][k]
            print
sys.exit()
