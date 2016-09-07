#!/usr/bin/env python

# electronicdos.py v0.7 03-06-2014 Jeff Doak jeff.w.doak@gmail.com
# Changelog:
# v0.6 - Added method to write multiple densities of states to the same file.
#          This method is untested as of 11/5/2013
# v0.5 - Added various l-, m-, and s-summation methods.
#      - Added non-collinear calculation flag. (Methods do not currently handle
#          non-collinear calculations, however.)
# v0.4 - Cleaned up sum_site_dos method.
# v0.3 - Adding to method sum_site_dos to include summing over several sets of
#        atoms concurrently without needing to call the method multiple times.
# v0.1 - Class name is the same as a class I use to integrate a DOS and
#        calculate finite temperature properties. I should either combine them
#        at some point, or rename one of the classes.
import numpy as np
import sys, subprocess, re
import aces.UnitCell.unitcell as uc

class ElectronicDOS:
    """
    Class to store and perform operations on an electronic density of states,
    as calculated by VASP. ElectronicDOS reads in a VASP-formatted DOSCAR file,
    including site- and lm-decomposed densities of states. The resulting
    densities of states can be summed over different sites and orbitals.

    Instance attributes of ElectronicDOS:
    - dosfile - name of DOSCAR-like file. Defaults to DOSCAR.
    - outfile - name of OUTCAR-like file. Defaults to OUTCAR.
    - posfile - name of POSCAR-like file. Defaults to POSCAR.
    - ispin   - flag to specify if DOS is spin-polarized (0 - not SP, 1 - SP).
    - lorbit  - flag to specify LORBIT setting of DOS file.
    - rwigs   - flag to specify whether RWIGS was set in calculation.
    - e_min   - minimum energy in DOSCAR.
    - e_max   - maximum energy in DOSCAR.
    - efermi  - Fermi energy as calculated by VASP.
    - n_atoms - number of atoms in unit cell of calculation.
    - n_dos   - number of energy points in electronic DOS.
    - energy  - numpy array containing list of energies at which DOS is 
                evaluated.
    - tot_dos - numpy array containing total electronic DOS. First column of
                array contains total DOS if not spin-polarized or spin-up DOS
                if spin-polarized. Second column of array contains zeros if not
                spin-polarized or -1*spin-down DOS if spin-polarized.
    - site_dos - numpy array containing atomic site- and lm-projected DOS.
                   Array is 3-dimensional. First index corresponds to the atom 
                   number, second index corresponds to energy of DOS, and third 
                   index corresponds to l, m, and spin DOS. If not
                   spin-polarized, last index goes in order l-m,...,l,...,l+m.
                   If spin-polarized, last index goes: l-m+up,l-m-dn,...,l+up,
                   l-dn,.../l+m+up,l+m-dn.
    """

    def __init__(self,dosfile="DOSCAR",outfile="OUTCAR",posfile="POSCAR"):
        self.dosfile = dosfile
        self.outfile = outfile
        self.posfile = posfile
        self.read_outcar()
        self.read_doscar()
        self.unit_cell = uc.UnitCell(self.posfile)

    def read_outcar(self):
        """
        Reads DOS parameters from OUTCAR file.
        """
        # Default values of OUTCAR parameters:
        self.ispin = 0
        self.lorbit = 0
        self.rwigs = 0
        self.noncoll = 0
        try:
            outcar = open(self.outfile,'r')
        except IOError:
            print "Error reading OUTCAR file."
            print "Default values for input settings will be used."
            return
        file = outcar.read()
        outcar.close()
        # Set ispin flag
        spin_str = r"ISPIN *= *([0-9]*)"
        spin_reg = re.compile(spin_str)
        if spin_reg.findall(file)[0] == '1':
            self.ispin = 0
        else:
            self.ispin = 1
        # Set lorbit flag
        orbit_str = r"LORBIT *= *([0-9]*)"
        orbit_reg = re.compile(orbit_str)
        self.lorbit = int(orbit_reg.findall(file)[0])
        # Set rwigs flag
        rwigs_str = r"RWIGS *= *[-0-9 .]*"
        rwigs_reg = re.compile(rwigs_str)
        temp = rwigs_reg.findall(file)[-1].split()
        temp.pop(0); temp.pop(0)
        for i in range(len(temp)):
            if i != '-1.00':   # Default RWIGS is a string of -1.00
                self.rwigs = 1
                break
        # Set noncoll flag
        noncoll_str = r"LNONCOLLINEAR *= *([FT])"
        noncoll_reg = re.compile(noncoll_str)
        if noncoll_reg.findall(file)[0] == "F":
            self.noncoll = 0
        else:
            self.noncoll = 1
        return

    def read_doscar(self):
        """
        Reads in a doscar file to grab the density of states as a function of
        energy.
        """
        try:
            input_ = open(self.dosfile,'r')
        except IOError:
            print "Error reading "+dosfile+" file."
            print "Program will now exit."
            sys.exit(1)
        # Read then discard header information
        self.n_atoms = int(input_.readline().split()[0])
        for i in range(4):
            input_.readline()
        # Read in Fermi Energy
        line = input_.readline().split()
        self.e_max = float(line[0])
        self.e_min = float(line[1])
        self.n_dos = int(line[2])
        self.e_fermi = float(line[3])
        # Read in total electronic DOS
        energy = []; tot_dos = []
        for i in range(self.n_dos):
            tot_dos.append([])
            line = input_.readline().split()
            energy.append(float(line[0]))
            if self.ispin == 0:  # non-spin polarized or non-collinear calc
                tot_dos[i].append(float(line[1]))  # DOS includes spin up and down
                #tot_dos[i].append(0.0)
            else:  # spin-polarized calculation
                tot_dos[i].append(float(line[1]))
                tot_dos[i].append(-1.*float(line[2]))
        self.energy = np.array(energy)
        #self.tot_dos = np.array(tot_dos)/float(self.n_atoms)
        self.tot_dos = np.array(tot_dos)
        # Read in site-projected electronic DOS.
        if (self.lorbit >= 10) or (self.rwigs == 1 and self.lorbit < 10):
            site_dos = []
            # Loop over each atom in the calculation
            for i in range(self.n_atoms):
                site_dos.append([])
                input_.readline()  # Discard site-projected header line
                for j in range(len(self.energy)):
                    site_dos[i].append([])
                    line = input_.readline().split()
                    for k in range(1,len(line)):
                        if self.ispin == 0:
                            site_dos[i][j].append(float(line[k]))
                        else:
                            site_dos[i][j].append((-1.)**(k-1)*float(line[k]))
            self.site_dos = np.array(site_dos)

    def sum_site_dos(self,list_,dos=None):
        """
        Sums site-projected DOS over a list of atoms. If list_ is a 2-d list,
        each set of atoms will be summed over seperately, and the returned array
        will have the same dimensionality as the input dos. If list_ is a 1-d
        list, the returned array will have one less dimension than the input
        dos.
        """
        # Use site- and orbital-decomposed DOS if no input DOS given
        if dos == None:
            dos = self.site_dos
        # Determine shape of dos, and set shape of summed_dos to be 1 dimension
        # less.
        #print "Shape of list_ is: "+str(np.shape(list_))
        #print "Shape of dos is: "+str(np.shape(dos))
        if len(np.shape(dos)) == 3:
            #print "dos is 3-d array"
            a = len(dos[0,0])
            summed_dos = np.zeros((len(list_),self.n_dos,a))
        else:
            #print "dos is 2-d array"
            summed_dos = np.zeros((len(list_),self.n_dos))
        #print "Shape of summed_dos is: "+str(np.shape(summed_dos))
        # Assume list_ is a 1-d array 
        #print "Shape of summed_dos[0] is: "+str(np.shape(summed_dos[0]))
        #print "Shape of dos[0] is: "+str(np.shape(dos[0]))
        for i in range(len(list_)):
            try:
                len(list_[i])
                # list_[i] is an array
                for j in list_[i]:
                    summed_dos[i] += dos[j]
            except TypeError:
                # list_[i] is not an array
                array_flag = False
                summed_dos[i] = dos[list_[i]]
        return summed_dos

    def sum_lms_dos(self):
        """
        Sums l-, m-, and s-projected DOS for each site-projected DOS.
        Returns an array cotaining the total DOS proejcted onto each atom.
        Should work for collinear non-spin-polarized, collinear spin-polarized,
        and non-collinear calculations.
        """
        summed_dos = np.zeros((self.n_atoms,self.n_dos))
        if self.noncoll == 0:
            if self.ispin == 0:
                # collinear, non-spin-polarized
                for i in range(self.n_atoms):
                    for j in range(self.n_dos):
                        for k in range(len(self.site_dos[0,0])):
                            summed_dos[i,j] += self.site_dos[i,j,k]
            else:
                # collinear, spin-polarized
                for i in range(self.n_atoms):
                    for j in range(self.n_dos):
                        for k in range(len(self.site_dos[0,0])):
                            summed_dos[i,j] += self.site_dos[i,j,k]
        else:
            # non-collinear
            for i in range(self.n_atoms):
                for j in range(self.n_dos):
                    for k in range(0,len(self.site_dos[0,0]),4):
                        summed_dos[i,j] += self.site_dos[i,j,k]
        return summed_dos

    def sum_lm_dos(self):
        """
        Sums l- and m-projected DOS for each site-projected DOS, leaving the
        spin-up and spin-down differences. For collinear non-spin-polarized 
        calculations, half the DOS is plotted spin-up and half is plotted 
        spin-down. For non-collinear calculations, the total magnitiziation 
        density along the x-, y-, and z-axes are returned, with no change in
        sign (the sign of the returned magnitization density along x, y, and z
        are all positive).
        """
        if self.noncoll == 0:
            summed_dos = np.zeros((self.n_atoms,self.n_dos,2))
            if self.ispin == 0:
                # collinear, non-spin-polarized
                for i in range(self.n_atoms):
                    for j in range(self.n_dos):
                        for k in range(len(self.site_dos[0,0])):
                            summed_dos[i,j,0] += self.site_dos[i,j,k]/2.
                            summed_dos[i,j,1] -= self.site_dos[i,j,k]/2.
            else:
                # collinear, spin-polarized
                for i in range(self.n_atoms):
                    for j in range(self.n_dos):
                        for k in range(len(self.site_dos[0,0])):
                            if k%2 == 0:  # k is even, corresponding to spin-up DOS
                                summed_dos[i,j,0] += self.site_dos[i,j,k]
                            else:  # k is odd, corresponding to spin-down DOS
                                summed_dos[i,j,1] -= self.site_dos[i,j,k]
        else:
            # non-collinear 
            x = range(len(self.site_dos[0,0]))
            y = range(0,len(self.site_dos[0,0]),4)
            set = [ i for i in x if i not in y ]
            summed_dos = np.zeros((self.n_atoms,self.n_dos,3))
            for i in range(self.n_atoms):
                for j in range(self.n_dos):
                    for k in set:
                        summed_dos[i,j,k%4-1] += self.site_dos[i,j,k]
        return summed_dos

    def determine_l_list(self):
        """
        Determines now many l-orbitals are in the site-projected DOS.
        This method assumes that the # of l-orbitals is the same for each atom
        in the DOS.
        """
        if self.noncoll == 0:
            if self.ispin == 0:
                # collinear, non-spin-polarized
                if len(self.site_dos[0,0]) == 1:  # s-orbitals only
                    l_list = [1]
                elif len(self.site_dos[0,0]) == 4:  # s- and p-orbitals
                    l_list = [1,3]
                elif len(self.site_dos[0,0]) == 9:  # s-, p-, and d-orbitals
                    l_list = [1,3,5]
                elif len(self.site_dos[0,0]) == 16: # s-, p-, d-, and f-orbitals
                    l_list = [1,3,5,7]
                else:
                    print "Unexpected number of lm-orbitals found."
                    print "Program will now quit."
                    sys.exit(1)
            else:
                # collinear, spin-polarized
                if len(self.site_dos[0,0]) == 2:  # s-orbitals only
                    l_list = [2]
                elif len(self.site_dos[0,0]) == 8:  # s- and p-orbitals
                    l_list = [2,6]
                elif len(self.site_dos[0,0]) == 18:  # s-, p-, and d-orbitals
                    l_list = [2,6,10]
                elif len(self.site_dos[0,0]) == 32:  # s-, p-, d-, and f-orbitals
                    l_list = [2,6,10,14]
                else:
                    print "Unexpected number of lms-orbitals found."
                    print "Program will now quit."
                    sys.exit(1)
        else:
            # non-collinear
            if len(self.site_dos[0,0]) == 4:  # s-orbitals only
                l_list = [4]
            elif len(self.site_dos[0,0]) == 16:  # s- and p-orbitals
                l_list = [4,12]
            elif len(self.site_dos[0,0]) == 36:  # s-, p-, and d-orbitals
                l_list = [4,12,20]
            elif len(self.site_dos[0,0]) == 64:  # s-, p-, d-, and f-orbitals
                l_list = [4,12,20,28]
            else:
                print "Unexpected number of lms-orbitals found."
                print "Program will now quit."
                sys.exit(1)
        return l_list

    def sum_ms_dos(self):
        """
        Sums m- and s-projected DOS for each site-projected DOS.
        Returns an array containing the atom-, l-projected DOS.
        Works for non-collinear calculations.
        """
        l_list = self.determine_l_list()
        # Perform the summation
        summed_dos = np.zeros((self.n_atoms,self.n_dos,len(l_list)))
        if self.noncoll == 0:
            if self.ispin == 0:
                # collinear non-spin-polarized
                for i in range(self.n_atoms):
                    for j in range(self.n_dos):
                        n = 0
                        for l in range(len(l_list)):
                            for m in range(l_list[l]):
                                summed_dos[i,j,l] += self.site_dos[i,j,n]
                                n += 1
            else:
                # collinear spin-polarized
                for i in range(self.n_atoms):
                    for j in range(self.n_dos):
                        n = 0
                        for l in range(len(l_list)):
                            for m in range(l_list[l]):
                                summed_dos[i,j,l] += self.site_dos[i,j,n]
                                n += 1
        else:
            # non-collinear
            for i in range(self.n_atoms):
                for j in range(self.n_dos):
                    for l in range(len(l_list)):
                        for m in range(0,l_list[l],4):
                            summed_dos[i,j,l] += self.site_dos[i,j,4*l+m]
        return summed_dos

    def sum_m_dos(self):
        """
        Sums the m-projected DOS for each site-projected DOS. For
        non-spin-polarized calculations, the spin-up and spin-down contributions
        to each l-orbit are 1/2 total density of each l-orbit. For non-collinear
        calculations, the atom- and l-projected magnitization density along the
        x-, y-, and z-axes are summed over each m quantum number.
        """
        l_list = self.determine_l_list()
        # Perform summation
        if self.noncoll == 0:
            summed_dos = np.zeros((self.n_atoms,self.n_dos,2*len(l_list)))
            if self.ispin == 0:
                # collinear non-spin-polarized
                for i in range(self.n_atoms):
                    for j in range(self.n_dos):
                        n = 0
                        for l in range(len(l_list)):
                            for m in range(l_list[l]):
                                summed_dos[i,j,2*l] += self.site_dos[i,j,n]/2.
                                summed_dos[i,j,2*l+1] -= self.site_dos[i,j,n]/2.
                                n += 1
            else:
                # collinear spin-polarized
                for i in range(self.n_atoms):
                    for j in range(self.n_dos):
                        n = 0
                        for l in range(len(l_list)):
                            for m in range(l_list[l]):
                                summed_dos[i,j,2*l+n%2] += self.site_dos[i,j,n]*(-1.)**(n)
                                n += 1
        else:
            # non-collinear
            summed_dos = np.zeros((self.n_atoms,self.n_dos,3*len(l_list)))
            for i in range(self.n_atoms):
                for j in range(self.n_dos):
                    for l in range(len(l_list)):
                        x = range(l_list[l])
                        y = range(0,l_list[l],4)
                        set = [ i for i in x if i not in y ]
                        for m in set:
                            summed_dos[i,j,3*l+m%3-1] += self.site_dos[i,j,4*l+m]
        return summed_dos
        
    def shift_energy(self,energy):
        """
        Adds the constant energy to the energy scale.
        """
        self.e_min += energy
        self.e_max += energy
        self.e_fermi += energy
        self.energy += energy

    def scale_dos(self,scale,dos=None):
        """
        Returns dos scaled by the factor scale.
        """
        if dos == None:
            dos = self.tot_dos
        return dos*scale

    def get_band_gaps(self,dos=None,tol=1e-3):
        """
        Returns a list of tuples containing the lower and upper bounds of each
        gap in the supplied dos. If no dos supplied, the total dos is used.
        """
        if dos == None:
            dos = self.tot_dos
        else:
            dos = np.array(dos)
        gaps = []
        if np.shape(dos)[0] == self.n_dos:
            # First entry of dos is list of energy points
            if len(np.shape(dos)) == 3:
                for i in range(len(dos[0])):
                    gaps.append([])
                    for j in range(len(dos[0,i])):
                        gaps.append([])
                        # Look at first energy of DOS 
                        if np.abs(dos[0,i,j]) < tol:
                            # DOS starts in a gap
                            flag = 1
                            #gaps[i][j].append([self.energy[0]])
                            gaps[i].append([])
                        else:
                            # DOS starts in a band
                            flag = 0
                        # Loop over all DOS points except first and last
                        for k in range(1,self.n_dos):
                            if (np.abs(dos[k,i,j]) > tol) and (flag == 1):
                                # Found upper bound of a gap
                                flag = 0
                                gaps[i][j][-1].append(self.energy[k])
                            elif (np.abs(dos[k,i,j]) < tol) and (flag == 0):
                                # Found lower bound of a gap
                                flag = 1
                                gaps[i][j].append([self.energy[k-1]])
            elif len(np.shape(dos)) == 2:
                for i in range(len(dos[0])):
                    gaps.append([])
                    if np.abs(dos[0,i]) < tol:
                        # DOS starts in a gap
                        flag = 1
                        #gaps[i].append([self.energy[0]])
                        gaps[i].append([])
                    else:
                        # DOS starts in a band
                        flag = 0
                    for k in range(1,self.n_dos):
                        if (np.abs(dos[k,i]) > tol) and (flag == 1):
                            # Found upper bound of a gap
                            flag = 0
                            gaps[i][-1].append(self.energy[k])
                        elif (np.abs(dos[k,i]) < tol) and (flag == 0):
                            # Found lower bound of a gap
                            flag = 1
                            gaps[i].append([self.energy[k-1]])
            else:
                print "Shape of dos is unexpected."
                print "Program will now quit!"
                sys.exit(1)
        elif np.shape(dos)[1] == self.n_dos:
            # First entry of dos is list of sites
            # Second entry of dos is list of energy points
            if len(np.shape(dos)) == 3:
                for i in range(len(dos)):
                    gaps.append([])
                    for j in range(len(dos[i,0])):
                        gaps.append([])
                        if np.abs(dos[i,0,j]) < tol:
                            # DOS starts in a gap
                            flag = 1
                            #gaps[i][j].append([self.energy[0]])
                            gaps[i].append([])
                        else:
                            # DOS starts in a band
                            flag = 0
                        for k in range(1,self.n_dos):
                            if (np.abs(dos[i,k,j]) > tol) and (flag == 1):
                                # Found upper bound of a gap
                                flag = 0
                                gaps[i][j][-1].append(self.energy[k])
                            elif (np.abs(dos[i,k,j]) < tol) and (flag == 0):
                                # Found lower bound of a gap
                                flag = 1
                                gaps[i][j].append([self.energy[k-1]])
            elif len(np.shape(dos)) == 2:
                for i in range(len(dos)):
                    gaps.append([])
                    if np.abs(dos[i,0]) < tol:
                        # DOS starts in a gap
                        flag = 1
                        #gaps[i].append([self.energy[0]])
                        gaps[i].append([])
                    else:
                        # DOS starts in a band
                        flag = 0
                    for k in range(1,self.n_dos):
                        if (np.abs(dos[i,k]) > tol) and (flag == 1):
                            # Found upper bound of a gap
                            flag = 0
                            gaps[i][-1].append(self.energy[k])
                        elif (np.abs(dos[i,k]) < tol) and (flag == 0):
                            # Found lower bound of a gap
                            flag = 1
                            gaps[i].append([self.energy[k-1]])
            else:
                print "Shape of dos is unexpected."
                print "Program will now quit!"
                sys.exit(1)
        else:
            # dos has a strange/incorrect formatting!
            print "Unexpected formatting for dos file."
            print "Program will now quit!"
            sys.exit(1)
        return np.array(gaps)

    def write_dos(self,dos=None):
        """
        Returns a string containing an electronic DOS or DOS's formatted for
        plotting.
        """
        if dos is None:
            dos = self.site_dos
        else:
            dos = np.array(dos)
        output = ""
        if np.shape(dos)[0] == self.n_dos:
            # First entry of dos is list of energy points
            for i in range(self.n_dos):
                output += str(self.energy[i])
                for j in range(len(dos[i])):
                    output += " "+str(dos[i,j])
                output += "\n"
        elif np.shape(dos)[1] == self.n_dos:
            # First entry of dos is list of sites
            if len(np.shape(dos)) == 3:
                for i in range(self.n_dos):
                    output += str(self.energy[i])
                    for j in range(len(dos)):
                        for k in range(len(dos[j,i])):
                            output += " "+str(dos[j,i,k])
                    output += "\n"
            elif len(np.shape(dos)) == 2:
                for i in range(self.n_dos):
                    output += str(self.energy[i])
                    for j in range(len(dos)):
                        output += " "+str(dos[j,i])
                    output += "\n"
        else:
            # dos has a strange/incorrect formatting!
            print "Unexpected formatting for dos file."
            print "Program will now quit!"
            sys.exit(1)
        return output

    def write_doses(self,dos_list):
        """
        Returns a string containing a set of electronic DOS's formatted for
        plotting.
        """
        output = ""
        for e in range(self.n_dos):
            output += str(self.energy[e])
            for i in range(len(dos_list)):
                if np.shape(dos_list[i])[0] == self.n_dos:  
                    # First energy of dos_list[i] is list of energy points
                    for j in range(len(dos_list[i][e])):
                        output += " "+str(dos_list[i][e,j])
                elif np.shape(dos_list[i])[1] == self.n_dos:
                    # First entry of dos_list[i] is list of sites
                    if len(np.shape(dos_list[i])) == 3:
                        for j in range(len(dos_list[i])):
                                for k in range(len(dos_list[i][j,e])):
                                    output += " "+str(dos_list[i][j,e,k])
                    elif len(np.shape(dos_list[i])) == 2:
                        for j in range(len(dos_list[i])):
                            output += " "+str(dos_list[i][j,e])
            output += "\n"
        return output

    def generate_z_list(self,tol=1e-8):
        """
        Returns a 2-d list of each plane of atoms, where the first index of list
        contains each plane of atoms along the z-axis, and the second index
        contains a list of atoms in each of the planes.
        """
        z_positions = self.unit_cell.atom_positions[:,2]
        z_sorted = np.argsort(z_positions)
        result = []
        result.append([z_sorted[0]])
        for i in range(1,len(z_sorted)):
            if abs(z_positions[z_sorted[i]] - z_positions[z_sorted[i-1]]) < tol:
                result[-1].append(z_sorted[i])
            else:
                result.append([z_sorted[i]])
        return np.array(result)

    def make_z_dos(self):
        """
        Sums site-projected DOS over lists of atoms lying in the same plane
        along the z-axis of the unit cell.
        """
        summed_dos = self.sum_lms_dos()
        z_lists = self.generate_z_list()
        z_dos = []
        for i in range(len(z_lists)):
            z_dos.append(self.sum_site_dos(z_lists[i],summed_dos))
        return np.array(z_dos)

    def make_half_dos(self):
        """
        Sums site-projected DOS over two halves of the unit cell divided in two
        along the z-direction.
        """
        summed_dos = self.sum_lms_dos()
        z_lists = self.generate_z_list()
        half = len(z_lists)/2
        half_dos = []
        half_dos.append(self.sum_site_dos(z_lists[0:half].flatten(),summed_dos))
        half_dos.append(self.sum_site_dos(z_lists[half:].flatten(),summed_dos))
        return np.array(half_dos)
    
    def dos_spline(self,dos=None):
        """
        Fits a DOS to a cubic spline. Useful for adding/subtracting DOS curves
        that have not been sampled on the same set of energies.
        """
        from scipy.interpolate import InterpolatedUnivariateSpline
        if dos == None:
            dos = self.tot_dos[:,0]
        splinedos = InterpolatedUnivariateSpline(self.energy,dos)
        return splinedos

    def dos_difference(self,dos1,dos2):
        """
        Subtracts two DOS curves. Assumes that dos1 and dos2 are cubic-splines.
        """
        difference = dos1(self.energy) - dos2(self.energy)
        return difference

if __name__ == "__main__":
    dos = ElectronicDOS()
    dos.shift_energy(-1.*dos.e_fermi)
    print np.shape(dos.tot_dos)
    #print dos for each half of the interface
    #sorted = dos.make_half_dos()
    #sorted = dos.scale_dos(1./dos.n_atoms,sorted)
    #dos_text = dos.write_dos(sorted)
    sorted = dos.make_z_dos()
    gaps = dos.get_band_gaps(dos=sorted,tol=1e-5)
    print np.shape(gaps)
    for i in range(len(gaps)):
        for j in range(len(gaps[i])):
            for k in range(len(gaps[i][j])):
                print i,gaps[i][j][k]
                print str(i+1),gaps[i][j][k]
                print

