#!/usr/bin/env python

# d_band_alignment.py v1.0 11-06-2013 Jeff Doak jeff.w.doak@gmail.com

from electronicdos import *
import sys


defdir = str(sys.argv[1])
purdir = str(sys.argv[2])

def_dos = ElectronicDOS(defdir+"/DOSCAR",defdir+"/OUTCAR",defdir+"/POSCAR")
pur_dos = ElectronicDOS(purdir+"/DOSCAR",purdir+"/OUTCAR",purdir+"/POSCAR")

tol = 1e-9
for i in range(def_dos.n_dos):
    def_energy = def_dos.energy[i]
    def_val = def_dos.tot_dos[i]
    if def_val > tol:
        break
for i in range(pur_dos.n_dos):
    pur_energy = pur_dos.energy[i]
    pur_val = pur_dos.tot_dos[i]
    if pur_val > tol:
        break

print "$E_{d-band}^{pure} (eV)$","$E_{d-band}^{def} (eV)$","$\Delta E_{d-band} (eV)$"
print pur_energy,def_energy,def_energy - pur_energy
sys.exit()
