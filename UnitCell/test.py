#!/usr/bin/env python
import re
import sys
import numpy as np

# Regex for gulp unit cell vectors at start of output file
reg_top = []
reg_a_top = re.compile(r"\ba *= *([0-9]*[.][0-9]*)")
reg_b_top = re.compile(r"\bb *= *([0-9]*[.][0-9]*)")
reg_c_top = re.compile(r"\bc *= *([0-9]*[.][0-9]*)")
reg_alpha_top = re.compile(r"alpha *= *([0-9]*[.][0-9]*)")
reg_beta_top = re.compile(r"beta *= *([0-9]*[.][0-9]*)")
reg_gamma_top = re.compile(r"gamma *= *([0-9]*[.][0-9]*)")

# Regex for gulp unit cell vectors at end of output file
reg_a_bot = re.compile(r"\ba *([0-9]*[.][0-9]*) *Angstrom")
reg_b_bot = re.compile(r"\bb *([0-9]*[.][0-9]*) *Angstrom")
reg_c_bot = re.compile(r"\bc *([0-9]*[.][0-9]*) *Angstrom")
reg_alpha_bot = re.compile(r"alpha *([0-9]*[.][0-9]*) *Degrees")
reg_beta_bot = re.compile(r"beta *([0-9]*[.][0-9]*) *Degrees")
reg_gamma_bot = re.compile(r"gamma *([0-9]*[.][0-9]*) *Degrees")

# Regex for checking if optimisation was performed
reg_opti = re.compile(r"optimise *- *perform optimisation run")
reg_num_atoms = re.compile(r"Total number atoms/shells *= *([0-9]*)")

# Regex for atomic positions at bottom of output file
reg_start = re.compile(r"Final fractional coordinates of atoms")
reg_start2 = re.compile(r"Label *\(Frac\) *\(Frac\) *\(Frac\) *\(Angs\) *\n-*")
reg_stop2 = re.compile(r"-*\n\n *Final Cartesian lattice vectors")
reg_stop = re.compile(r"Final Cartesian lattice vectors")

# Regex for atomic positions at top of output file
reg_start_top = re.compile(
                r"Label *\(Frac\) *\(Frac\) *\(Frac\) *\(e\) *\(Frac\) *\n-*\n")
reg_stop_top = re.compile(r"-*\n*\**\n\* *General input information")

# Read in gulp output file
gulp_in = open("gulp.in","r")
gulp_in.readline()
cell_vec = np.zeros((3,3))
num_atoms = 0
atom_names = []
atom_positions = []
line = gulp_in.readline()
for line in gulp_in:
    if line.count("cell"):
        line = gulp_in.readline().split()
        a = float(line[0])
        b = float(line[1])
        c = float(line[2])
        alpha = float(line[3])
        beta = float(line[4])
        gamma = float(line[5])
        print a,b,c,alpha,beta,gamma
    elif line.count("title"):
        name = gulp_in.readline()
        gulp_in.readline()  # Discard line containing 'end'
    elif line.count("vector"):
        for i in range(3):
            line = gulp_in.readline().split()
            cell_vec[i,0] = line[0]
            cell_vec[i,1] = line[1]
            cell_vec[i,2] = line[2]
    elif line.count("frac"):
        while True:
            line = gulp_in.readline().split()
            if len(line) != 4:
                break




sys.exit()

# Test lat param grabbing

a = reg_a_bot.search(lines)
n = reg_num_atoms.search(lines)
print n


#b = lines[start:stop]
#print b
start2 = reg_start2.search(lines).end()
stop2 = reg_stop2.search(lines).start()
#print start2
#print stop2
start = reg_start_top.search(lines).end()
stop = reg_stop_top.search(lines).start()

atom_lines = lines[start2+1:stop2-1].split('\n')
atom_names = []
num_atoms = len(atom_lines)
atom_positions = np.zeros((num_atoms,3))
for i in range(len(atom_lines)):
    line = atom_lines[i].split()
    atom_names.append(line[1])
    if len(line) <= 8:
        atom_positions[i,0] = float(line[3])
        atom_positions[i,1] = float(line[4])
        atom_positions[i,2] = float(line[5])
    else:
        atom_positions[i,0] = float(line[3])
        atom_positions[i,1] = float(line[5])
        atom_positions[i,2] = float(line[7])
#print atom_positions
print atom_names
print set(atom_names)
num_atom_types = len(set(atom_names))
atom_types = []
for i in range(num_atom_types):
    atom_types.append(
            atom_names.count(list(set(atom_names))[i]))
print atom_types
