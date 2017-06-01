#!/usr/bin/env python

# plot_partial_orbital_edos.py v1.1 09-23-2013 Jeff Doak jeff.w.doak@gmail.com

# This script sums the site-projected electronic DOS over each type of atom,
# adding each orbital of each atom type separately. The result is an atom_type-
# and orbital-projected electronic DOS containing a 3-d array of atom type,
# energy level, and orbital.

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
# Sum the densities of states over each orbital for each atom.
orbital_dos = doscar.sum_ms_dos()
# Create a list of each atom type to sum over.
type_list = []
n = 0
for i in range(len(doscar.unit_cell.atom_types)):
    type_list.append([])
    for j in range(doscar.unit_cell.atom_types[i]):
        type_list[i].append(n)
        n += 1
# Sum dos over sets of atoms.
partial_dos = doscar.sum_site_dos(type_list,orbital_dos)
print doscar.write_doses([doscar.tot_dos,partial_dos])
sys.exit()
