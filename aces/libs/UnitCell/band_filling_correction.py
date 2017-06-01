#!/usr/bin/env python

# band_filling_correction.py v1.0 6-12-2013 Jeff Doak jeff.w.doak@gmail.com
# Calculates the band-filling correction for a defect

from unitcell import *
import re,sys
import numpy as np

# 1st argument is directory containing defect calculation
# 2nd argument is directory containing reference calculation

ref_out = str(sys.argv[2])+"OUTCAR"

# Read in OUTCAR of defect calculation
def_out = open(str(sys.argv[1])+"OUTCAR","r")
#lines = def_out.readlines()
lines = def_out.read()
def_out.close()
def_kpts = []
def_weights = []
def_energies = []
def_occupancies = []
kpt_reg = re.compile(r" *kpoint *[0-9] *: [0-9]*[.][0-9]* *[0-9]*[.][0-9]* *[0-9][.][0-9]*")
print kpt_reg.findall(lines)
exit()
for i in range(len(lines)):
    if lines[i].startswith("  (the norm of the test charge is"):
        j = 1
        k = 0
        while k < poscar.num_atoms:
            line = lines[i+j].split()
            while len(line) > 0:
                k = int(line.pop(0))
                charged_pots[k-1] = float(line.pop(0))
            j += 1
        break

# Read in OUTCAR of neutral calculation
outcar = open(neutralout,"r")
lines = outcar.readlines()
outcar.close()
neutral_pots = np.zeros(perfect.num_atoms)
for i in range(len(lines)):
    if lines[i].startswith("  (the norm of the test charge is"):
        j = 1
        k = 0
        while k < perfect.num_atoms:
            line = lines[i+j].split()
            while len(line) > 0:
                k = int(line.pop(0))
                neutral_pots[k-1] = float(line.pop(0))
            j += 1
        break


# Find half the shortest defect-image distance
r_min = np.linalg.norm(poscar.cell_vec[0])/2. #huge assumption that the cell is cubic
r_min_neut = np.linalg.norm(perfect.cell_vec[0])/2. #huge assumption that the cell is cubic

r_min = poscar.get_volume()**(1/3.)/2.
r_min_neut = perfect.get_volume()**(1/3.)/2.

print r_min,r_min_neut

# Average electrostatic potential of charged cell
charged_radii = [ np.linalg.norm(i-center) for i in poscar.atom_positions ]
sphere_charged_pots = []
for i in range(len(charged_radii)):
    if charged_radii[i] >= r_min:
        sphere_charged_pots.append(charged_pots[i])
chg_avg = np.mean(sphere_charged_pots)

# Average electrostatic potential of neutral cell
neutral_radii = [ np.linalg.norm(i-neut_center) for i in perfect.atom_positions ]
sphere_neutral_pots = []
for i in range(len(neutral_radii)):
    if neutral_radii[i] >= r_min_neut:
        sphere_neutral_pots.append(neutral_pots[i])
neutral_avg = np.mean(sphere_neutral_pots)

# Calculate average electrostatic potential difference
delta_V_pa = chg_avg - neutral_avg
print "Average over atoms outside radius",r_min,"A centered around defect"
print "$\Delta V_{PA} (eV)$ &","# atoms in average &","Std. Dev. (eV)"
print delta_V_pa,len(sphere_charged_pots),len(sphere_neutral_pots)

if len(charged_pots) == len(neutral_pots):
    delta_pot = charged_pots - neutral_pots
    radii = [ np.linalg.norm(i-center) for i in poscar.atom_positions ]
    sphere_pots = []
    for i in range(len(radii)):
        if radii[i] >= r_min:
            sphere_pots.append(delta_pot[i])
    print ""
    print "Average over atoms outside radius",r_min,"A centered around defect"
    print "$\Delta V_{el} (eV)$ &","# atoms in average &","Std. Dev. (eV)"
    print np.mean(sphere_pots),len(sphere_pots),np.std(sphere_pots)
exit()


# Find atomic radii of each atom away from defect atom


# Calculate electrostatic potential alignment
delta_pot = charged_pots - neutral_pots
radii = [ np.linalg.norm(i-center) for i in poscar.atom_positions ]
sorted_radii = np.sort(radii)[::-1]
sorted_delta_pot = [ delta_pot[i] for i in np.argsort(radii) ][::-1]

# Average over the last 2 sets of radii - should correspond to furthest cation
# and anion from defect.
uniq_pots,uniq_index =  np.unique(sorted_delta_pot,return_index=True)
last_2_avg = np.mean([ sorted_delta_pot[i] for i in np.sort(uniq_index)[0:2] ])

# Find half the shortest defect-image distance
r_min = np.linalg.norm(poscar.cell_vec[0])/2. #huge assumption that the cell is cubic

# average over all atoms lying outside sphere with radius r_min
sphere_pots = []
for i in range(len(radii)):
    if radii[i] >= r_min:
        sphere_pots.append(delta_pot[i])



# various averages
print "$\Delta V_{el} (eV)$","# atoms in average"
print np.mean(sorted_delta_pot),poscar.num_atoms
print np.mean(sorted_delta_pot[0:poscar.num_atoms/2]),poscar.num_atoms/2
print np.mean(sorted_delta_pot[0:poscar.num_atoms/4]),poscar.num_atoms/4
print np.mean(sorted_delta_pot[0:poscar.num_atoms/6]),poscar.num_atoms/6
print np.mean(sorted_delta_pot[0:poscar.num_atoms/8]),poscar.num_atoms/8
print np.mean(sorted_delta_pot[0:poscar.num_atoms/10]),poscar.num_atoms/10
print np.mean(sorted_delta_pot[0:poscar.num_atoms/12]),poscar.num_atoms/12
print np.mean(sorted_delta_pot[0:poscar.num_atoms/14]),poscar.num_atoms/14
print np.mean(sorted_delta_pot[0:poscar.num_atoms/16]),poscar.num_atoms/16
print np.mean(sorted_delta_pot[0:poscar.num_atoms/18]),poscar.num_atoms/18
print np.mean(sorted_delta_pot[0:poscar.num_atoms/20]),poscar.num_atoms/20
print np.mean(sorted_delta_pot[0:2]),2
print last_2_avg,"last cation and anion"

print "Average over atoms outside radius",r_min,"A centered around defect"
print "$\Delta V_{el} (eV)$ &","# atoms in average &","Std. Dev. (eV)"
print np.mean(sphere_pots),len(sphere_pots),np.std(sphere_pots)

sys.exit()
