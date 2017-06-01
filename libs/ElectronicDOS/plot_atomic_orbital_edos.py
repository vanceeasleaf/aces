#!/usr/bin/env python

# plot_atomic_orbital_dos.py v1.0 10-03-2012 Jeff Doak jeff.w.doak@gmail.com
# Plots the site- and lms-decomposed electronic DOS. The shape of the printed
# array is (NEDOS,N_atoms,9 (or 18 if spin-polarized))

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
print doscar.write_dos(doscar.site_dos)
sys.exit()
