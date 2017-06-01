#!/usr/bin/env python

# plot_partial_edos.py v1.0 12-12-2012 Jeff Doak jeff.w.doak@gmail.com

# This script sums the site- and orbital-projected electronic DOS over each
# orbital and each type of atom. The result is a partial electronic DOS
# containing one list of density vs electron energy level for each type of atom
# listed in the POSCAR file.

from electronicdos import *
import sys

# Read in electronic DOS from DOCAR file
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

# Shift the energy scale so that the Fermi Energy is at 0 eV.
#doscar.shift_energy(-1.*doscar.e_fermi)
# Create a list of each atom type to sum over.
type_list = []
n = 0
for i in range(len(doscar.unit_cell.atom_types)):
    type_list.append([])
    for j in range(doscar.unit_cell.atom_types[i]):
        type_list[i].append(n)
        n += 1
# Sum dos over sets of atoms.
atom_dos = doscar.sum_lms_dos()
partial_dos = doscar.sum_site_dos(type_list,atom_dos)
#print doscar.write_dos(partial_dos)
print doscar.write_doses([doscar.tot_dos,partial_dos])
sys.exit()
