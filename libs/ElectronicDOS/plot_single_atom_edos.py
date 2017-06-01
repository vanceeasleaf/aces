#!/usr/bin/env python

# plot_single_edos.py v1.0 01-07-2013 Jeff Doak jeff.w.doak@gmail.com

# This script plots the site-projected electronic DOS of a set of atoms,
# keeping the orbital-projection of each atom intact. The result is an atom-
# and orbital-projected partial-electronic DOS containing a 3-d array of atom,
# energy level, and orbital.

from electronicdos import *
import sys

# Read in electronic DOS from DOCAR file
#if len(sys.argv) == 2:
#    dosdir = str(sys.argv[1])
#    doscar = ElectronicDOS(dosdir+"/DOSCAR",dosdir+"/OUTCAR",dosdir+"/POSCAR")
#elif len(sys.argv) == 4:
#    dosfile = str(sys.argv[1])
#    outfile = str(sys.argv[2])
#    posfile = str(sys.argv[3])
#    doscar = ElectronicDOS(dosfile,outfile,posfile)
#else:
#    doscar = ElectronicDOS()
atom_list = [int(i) for i in sys.argv[1:]]
print atom_list
sys.exit()

# Shift the energy scale so that the Fermi Energy is at 0 eV.
#doscar.shift_energy(-1.*doscar.e_fermi)
# Sum the densities of states over each orbital for each atom.
orbital_dos = doscar.sum_lms_dos()
# Create a list of each atom type to sum over.
type_list = []
n = 0
for i in range(len(doscar.unit_cell.atom_types)):
    type_list.append([])
    for j in range(doscar.unit_cell.atom_types[i]):
        type_list[i].append(n)
        n += 1
# Sum dos over sets of atoms.
partial_dos = doscar.sum_site_dos(type_list)
print doscar.write_dos(partial_dos)
sys.exit()
