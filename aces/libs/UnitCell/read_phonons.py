#!/usr/bin/python

# v0.2 3-12-2012 Jeff Doak jeff.w.doak@gmail.com
import numpy as np
import sys
from unitcell import *

def read_normal_modes(in_name="phonons.out"):
    """Reads in phonon normal mode coordinates from a phonons.out file."""
    in_file = open(in_name,"r")
    line = in_file.readline().split()
    num_atoms = int(line[0])
    num_q_points = int(line[1])
    in_file.readline()  # Discard q-point line
    num_modes = 3*num_atoms
    num_lines = num_modes/6  # Six frequencies per line
    freqs = []
    for i in range(num_lines):
        line = in_file.readline().split()
        for j in line:
            freqs.append(float(j))
    freqs = np.array(freqs)
    normal_modes = []
    for i in range(num_modes):
        normal_modes.append([])
        in_file.readline()  # Discard blank line
        for j in range(num_atoms):
            line = in_file.readline().split()
            for k in range(3):
                #normal_modes[i].append(complex(line[k+1],line[k+4]))
                normal_modes[i].append(float(line[k+1]))
    in_file.close()
    normal_modes = np.array(normal_modes)
    return freqs,normal_modes

def read_masses(mass_name="apos.dat"):
    """Reads in the atomic masses from the apos.dat file."""
    mass_file = open(mass_name,"r")
    # Discard unit cell parameters and supercell size.
    for i in range(4):
        mass_file.readline()
    # Get number of atom types
    num_atom_types = int(mass_file.readline().split()[0])
    # Get list of atomic masses (strings).
    masses = mass_file.readline().split()
    # Get list containing number of atoms of each type (strings).
    atom_types = mass_file.readline().split()
    mass_file.close()
    mass_vec = []
    for i in range(num_atom_types):
        for j in range(int(atom_types[i])):
            for k in range(3):
                mass_vec.append(float(masses[i]))
    mass_vec = np.array(mass_vec)
    return mass_vec #  numpy array with size #atomsx1.

def calc_atomic_displacements(normal_modes,mass_vec):
    """Calculates atomic displacements from phonon normal modes and the atomic
    masses."""
    normal_disp = normal_modes/np.sqrt(mass_vec)
    #normal_disp = normal_disp/np.linalg.norm(normal_disp)
    normal_disp = np.round(normal_disp,decimals=6)
    return normal_disp

def decompose_displacements(normal_modes,disp):
    """Calculates how much of a set of atomic displacements lies along each
    normal mode vector of the undisplaced cell."""
    # I need to solve the equation x*A=b => x=b*A**-1
    # Where A is a matrix and x,b are vectors.
    # However, there could be a problem, depending on how exactly I defined my
    # matrices. I need to make sure that every thing works out ok.
    #weights = np.dot(np.linalg.inv(normal_modes).transpose(),disp)
    weights = np.dot(disp,np.linalg.inv(normal_modes))
    return weights

def decompose_displacements2(normal_modes,disp,masses):
    weights = 0
    for i in range(len(disp)):
        pass



def main():
    file1 = open(str(sys.argv[1]),'r')
    file2 = open(str(sys.argv[2]),'r')
    cell_1 = UnitCell(file1)
    cell_2 = UnitCell(file2)
    disp = displacements(cell_1,cell_2,conv='C')
    disp = disp.flatten()
    if len(sys.argv) > 3:
        in_name = str(sys.argv[3])
    else:
        in_name = "phonons.out"
    if len(sys.argv) > 4:
        mass_name = str(sys.argv[4])
    else:
        mass_name = "apos.dat"
    freqs,normal_modes = read_normal_modes(in_name)
    mass_vec = read_masses(mass_name)
    normal_disp = calc_atomic_displacements(normal_modes,mass_vec)
    weights = decompose_displacements(normal_disp,disp)
    print np.array_repr(weights,precision=8,suppress_small=True)
    print np.linalg.norm(weights)
    #print np.array_repr(weights/np.linalg.norm(weights)*100,precision=4,suppress_small=True)
    #print np.array_str(weights/np.linalg.norm(weights)*100)

def test():
    mass_vec = read_masses()
    freqs,normal_modes = read_normal_modes()
    normal_disp = normal_modes/np.sqrt(mass_vec)
    inv_disp = np.linalg.inv(normal_disp)
    inv_modes = np.linalg.inv(normal_modes)*np.sqrt(mass_vec)
    trans_modes = normal_modes.transpose()
    print "inverse of normal modes:"
    print inv_modes
    print "transpose:"
    print trans_modes
    print "difference:"
    print inv_modes - trans_modes
    
    sys.exit()
    other_disp = np.zeros((len(normal_modes),len(normal_modes[0])))
    for i in range(len(normal_modes)):
        for j in range(len(normal_modes[i])):
            other_disp[i,j] = normal_modes[i,j]/np.sqrt(mass_vec[j])
    diff = normal_disp - other_disp
    print diff
    normal_disp = np.round(normal_disp,decimals=6)

def test2():
    file1 = open(str(sys.argv[1]),'r')
    file2 = open(str(sys.argv[2]),'r')
    cell_1 = UnitCell(file1)
    cell_2 = UnitCell(file2)
    disp = displacements(cell_1,cell_2,conv='C')




if __name__ == "__main__":    
    #main()
    test()
