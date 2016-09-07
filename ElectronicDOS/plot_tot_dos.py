#!/usr/bin/env python

# plot_tot_dos.py v1.0 9-24-2012 Jeff Doak jeff.w.doak@gmail.com

from electronicdos import *
import sys

if len(sys.argv) == 2:
    dosdir = str(sys.argv[1])
    doscar = ElectronicDOS(dosdir+"/DOSCAR",dosdir+"/OUTCAR",dosdir+"/POSCAR")
elif len(sys.argv) == 4:
    dosfile = str(sys.argv[1])
    outfile = str(sys.argv[2])
    posfile = str(sys.argv[3])
    doscar = ElectronicDOS(dosfile,outfile,posfile)
else:
    doscar = ElectronicDOS()

#doscar.shift_energy(-1.*doscar.e_fermi)
print doscar.write_dos(doscar.tot_dos)
sys.exit()
