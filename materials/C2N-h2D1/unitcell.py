#!/usr/bin/python

# UnitCell Class v2.0.0 Jeff Doak jeff.w.doak@gmail.com 09-26-2013
# Changelog:
# v2.0.0 - Added read_atat method.
# v1.9.8 - Added add_atom method. May be a bug in this method!
#        - Added delete_atom method.
#        - Added change_atom_type method.
# v1.9.7 - Fixed bugs in shift method.
#        - Fixed bugs in in_cell method.
#        - Added get_volume method.
# v1.9.6 - Removed CHGCAR methods. These now work with a separate class,
#          ChargeDensity
# v1.9.5 - Added Gulp input file input method. Currently only supports 'cell'
#          format, not 'vector' format.
# v1.9.3 - Added velocity flag to read_poscar method. With vel=True (default)
#          atom velocities will be read in. With vel=False, atom velocities will
#          be ignored. Useful for using read_poscar to read in the atomic
#          position portion of a CHGCAR-like file.
#        - Added read_chgcar method to read in a CHGCAR-like file. Assumes the
#          atomic positions have already been read in.
# v1.9.2 - Added method to perform a linear transformation on a unit cell. This
#          transformation alters the unit cell vectors, while leaving fractional
#          atomic positions unchanged.
#        - Fixed bug in read_poscar method where velocity attibutes were not
#          created if no velocities were present in the poscar file.
#        - Fixed bug in __init__ where self.atom_names was spelled incorrectly.
# v1.9.1 - Added (untested) section to __init__ to return a new UnitCell object
#          that is a copy of of a UnitCell object given as input_.
#        - Added instance method simple_supercell which takes three integers and
#          creates a supercell of the unit cell that is a*cell_vec[0],
#          b*cell_vec[1], and c*cell_vec[2] in size, and with the approriate
#          increase in number of atoms. This method is unfinshed!
# v1.9.0 - Added flag to displacement function to allow displacements to be
#          calculated in direct or cartesian coordinates. 
#        - Replaced list comprehension with re.search in the spacegroup method. 
#          Spacegroup method now removes the findsym.log file if findsym runs 
#          correctly.
#        - Reworked how conversion between direct and cartesian coordinates
#          happens. Made two package functions direct_to_cart and cart_to_direct,
#          which multiply a lattice and array together to convert the array from
#          direct to cart. coords. or vice versa. Changed the setter method for
#          self.convention to use these functions instead of the to_direct and
#          to_cart methods. 
#        - The read_poscar method now supports inputting MD velocities from the
#          POSCAR file. Created two new attributes, _vel_convention and
#          atom_velocities, which are analogous to the convention and
#          atom_positions attributes that currently exist. Added getter and setter
#          methods for the _vel_convention attribute.
#        - The output_vasp method now outputs MD velocities, if they are non-zero.
# v1.8   - Added method to read in gulp output file. Unit cell data can be read in
#          either from the final relaxed positions, or the initial positions.
#        - Modified output_ezvasp to only return a string, the user can deal with
#          appending the structure to an ezvasp file. 
#        - Added a function to convert a,b,c,alpha,beta,gamma to unit cell vectors
#          (not a class method, but a package function).
#        - Removed old commented out getters and setters.
#        - Commented out supercell method because it doesn't work as of now.
# v1.7   - Added method to move all atoms to their positions inside the unit 
#          cell, if they aren't already there. 
#        - Added a method to calculate the center of mass of a unit cell, reading
#          the masses from a POTCAR-formatted file.
#        - Fixed bug in set_convetion method where to_cart was named incorrectly.
# v1.6   - Added method to output ezvasp-formatted text or append to vasp.in file.
#        - Added get and set methods for changing the atom names. The variable
#          num_atom_types and atom_types can be broken by these methods still.
#        - Added method to change scale of unit cell in a way that does change the
#          lattice parameter and volume.
# v1.5   - Added method to calculate spacegroup using ISOTROPY package.
#        - Changed class method displacements to a function in unitcell package.
# v1.4   - Removed unneeded getters ans setters. Made remaining ones compatible
#          with python 2.4.3. 
#        - Edited displacement class method to return values instead of printing 
#          them, and added a second return format.
# v1.3   - Made class an extension of object; added getter and setter methods for
#          most attributes.
# v1.2   - Added shift method to shift atom positions by a constant vector.
# v1.1   - Added LAMMPS output method. Currently only supports orthogonal cells!

import numpy as np
import os,sys,subprocess,re

class UnitCell(object):
    """
    Class to read in and store crystal unit cell data.
    """

    # Define class attribute getters and setters.
    def get_scale(self):
        """
        Scale factor for unit cell vectors and atomic positions.
        """
        return self._scale

    def set_scale(self,new_scale=1.0):
        """
        Changes the scale of the unit cell to a new number, while leaving the
        unit cell properties untouched.
        """
        old_convention = self.convention
        self.convention = 'Direct'
        self.cell_vec = self.cell_vec*self._scale/float(new_scale)
        self._scale = new_scale
        self.set_convention(old_convention)

    def get_convention(self):
        """
        The convention used to store atomic positions. Only the first character
        in the string is used to determine convetion. C, c, K, or k indicate
        that atomic positions are stored in Cartesian coordiates. Anything else
        is stored in direct coordinates.
        """
        return self._convention

    def set_convention(self,string):
        cees = ['c','C','k','K']
        if not cees.count(string[0]):  # string was a form of Direct
            if cees.count(self._convention[0]):  # Atoms in cartesian
                #self.to_direct()
                self.atom_positions = cart_to_direct(
                        self.cell_vec,self.atom_positions)
                self._convention = 'Direct'
        else:  # string was a form of Cartesian
            if not cees.count(self._convention[0]):  # Atoms in direct
                #self.to_cart()
                self.atom_positions = direct_to_cart(
                        self.cell_vec,self.atom_positions)
                self._convention = 'Cartesian'

    def get_vel_convention(self):
        """
        The convention used to store atomic velocities. Only the first character
        in the string is used to determine convetion. C, c, K, or k indicate
        that atomic positions are stored in Cartesian coordiates. Anything else
        is stored in direct coordinates.
        """
        return self._vel_convention

    def set_vel_convention(self,string):
        vcees = ['c','C','k','K','\n','\t',' ']
        if not vcees.count(string[0]):  # string was Direct
            print "found direct velocities"
            print string
            if vcees.count(self._vel_convention[0]):  # Vel. in cartesian
                self.atom_velocities = cart_to_direct(
                        self.cell_vec,self.atom_velocities)
                self._vel_convention = 'Direct'
        else:  # string was Cartesian
            print "found cartesian velocities"
            print string
            if not vcees.count(self._vel_convention[0]):  # Vel. in direct
                self.atom_velocities = direct_to_cart(
                        self.cell_vec,self.atom_velocities)
                self._vel_convention = 'Cartesian'

    def get_atom_names(self):
        return self._atom_names

    def set_atom_names(self,list_):
        """
        Changes the names of atoms without allowing the length of the list
        containing the atom names to change.
        """
        if isinstance(list_,str):
            list_ = list_.split()
        if len(list_) <= self.num_atom_types:
            n = 0
            for i in range(len(list_)):
                for j in range(self.atom_types[i]):
                    self._atom_names[n] = list_[i]
                    n += 1
        # This method can change the number of atom types, without changing the
        # attribute self.num_atom_types!
        elif len(list_) == self.num_atoms:
            print "Warning: Number of atom types may have changed!"
            for i in range(len(list_)):
                self._atom_names[i] = list_[i]
        else:
            print "Warning: Length of names didn't match number of atom types!"
            n = 0
            for i in range(self.num_atom_types):
                for j in range(self.atom_types[i]):
                    self._atom_names[n] = list_[i]
                    n += 1

    scale = property(get_scale,set_scale)
    convention = property(get_convention,set_convention)
    vel_convention = property(get_vel_convention,set_vel_convention)
    atom_names = property(get_atom_names,set_atom_names)

    def __init__(self, input_=None,format_=None):
        if isinstance(input_,UnitCell):
            # Return a new UnitCell object that is a copy of input_
            self.name = input_.name
            self._scale = input_.scale
            self.cell_vec = np.array(input_.cell_vec)
            self.num_atom_types = input_.num_atom_types
            self.atom_types = input_.atom_types  # Will this return a copy or the reference?
            self._convention = input_.convention
            self.num_atoms = input_.num_atoms
            self.atom_positions = np.array(input_.atom_positions)
            self._atom_names = input_.atom_names  # Will this return a copy or the reference?
            self._vel_convention = input_.vel_convention
            self.atom_velocities = np.array(input_.atom_velocities)
        elif isinstance(input_,str):
            try:
                input_ = open(input_,'r')
            except IOError:
                print "Error reading input file."
                print "Empty UnitCell will be returned."
                input_ = None
        #if isinstance(input_,file):
        #else: # Assume input_ is a file
        elif input_ is not None:
            if format_ == None:
                self.read_poscar(input_)
            elif format_ == "lammps":
                self.read_lammps_dump(input_)
            elif format_ == "gulp_output":
                self.read_gulp_output(input_)
            elif format_ == "gulp_input":
                self.read_gulp_input(input_)
            elif format_ == "atat":
                self.read_atat(input_)
            elif format_ == "chgcar":
                self.read_poscar(input_,vel=False)
                self.read_chgcar(input_)
                self.scale_density()
            elif format_ == "locpot":
                self.read_poscar(input_,vel=False)
                self.read_chgcar(input_)
            else:
                print "Unknown format."
                print "Empty UnitCell will be returned."
                input_ = None
        if input_ == None:
            # Create 'empty' simple cubic unit cell as default.
            self.name = ""
            self._scale = 1.0
            self.cell_vec = np.identity(3)
            self.num_atom_types = 0
            self.atom_types = []
            self._atom_type_names = None
            self._convention = 'Direct'
            self.num_atoms = 0
            self.atom_positions = np.zeros((self.num_atoms,3))
            self._atom_names = [None]
            self._vel_convention = 'Cartesian'
            self.atom_velocities = np.zeros((self.num_atoms,3))

    def to_direct(self):
        """
        Convert atomic positions from cartesian coordinates to direct
        coordinates.
        """
        cees = ['c','C','k','K']
        if cees.count(self._convention[0]):
            for i in range(self.num_atoms):
                self.atom_positions[i] = np.dot(
                    np.linalg.inv(
                        self.cell_vec.transpose()),self.atom_positions[i])
            self._convention = 'Direct'

    def to_cart(self):
        """
        Convert atomic positions from direct coordinates to cartesian
        coordinates.
        """
        cees = ['c','C','k','K']
        if not cees.count(self._convention[0]):
            for i in range(self.num_atoms):
                self.atom_positions[i] = np.dot(
                    self.cell_vec.transpose(),self.atom_positions[i])
            self._convention = 'Cartesian'

    def shift(self,s_x,s_y,s_z,convention):
        """
        Shifts all atoms in the unit cell by s_i in the i-direction for
        x,y,z. convention determines whether shift is done to atoms in cartesian
        or direct coordinate systems.
        """
        shift = np.array((s_x,s_y,s_z))
        self.convention = convention
        for i in range(self.num_atoms):
            self.atom_positions[i] = self.atom_positions[i] + shift

    def change_scale(self,new_scale):
        """
        Changes the scale (i.e. lattice parameter and volume) of the unit cell
        to that of new_scale.
        """
        self._scale = float(new_scale)

    def in_cell(self):
        """
        Moves all atoms to their images within the unit cell.
        """
        self.set_convention("direct")
        for i in range(self.num_atoms):
            for j in range(3):
                if self.atom_positions[i,j] > 1.0:
                    self.atom_positions[i,j] = self.atom_positions[i,j] - 1.0
                elif self.atom_positions[i,j] < 0.0:
                    self.atom_positions[i,j] = 1.0 + self.atom_positions[i,j]

    def get_volume(self):
        """
        Returns the volume of the unitcell.
        """
        temp = self.scale
        self.scale = 1.0
        vol = np.dot(self.cell_vec[0],np.cross(self.cell_vec[1],self.cell_vec[2]))
        self.scale = temp
        return vol

    def add_atom(self,pos,type,convention,vel=np.zeros(3)):
        """
        Adds an atom of a given type to the unitcell at a position specified in
        either direct or cartesian coordinates. Optionally, a velocity for the
        new atom may be specified. The new atom will be added to the end of the
        list of atoms.
        """

        def uniq(list_):
            """
            Function to get all unique entries in a list, keeping entry order.
            """
            seen = set()
            seen_add = seen.add
            return [ x for x in list_ if x not in seen and not seen_add(x) ]

        self.convention = convention
        self.num_atoms += 1
        pos = np.array([pos])
        vel = np.array([vel])
        try:
            # check if type is an integer
            int(type)
        except ValueError:
            # type is a string
            name = type
            if name in self.atom_names:
                name_list = uniq(self.atom_names)
                for i in range(len(name_list)):
                    if name_list[i] == name:
                        type = i
                loc = 0
                for i in range(0,type+1):
                    loc += self.atom_types[i]
                self.atom_positions = np.insert(self.atom_positions,loc,pos,0)
                self.atom_velocities = np.insert(self.atom_velocities,loc,vel,0)
                self.atom_names.insert(loc,name)
                self.atom_types[type] += 1
            else:
                self.atom_positions = np.append(self.atom_positions,pos,0)
                self.atom_velocities = np.append(self.atom_velocities,vel,0)
                self.num_atom_types += 1
                self.atom_names.append(name)
                self.atom_types.append(1)
        else:
            # type is an integer
            if type < self.num_atom_types:
                loc = 0
                for i in range(0,type+1):
                    loc += self.atom_types[i]
                self.atom_positions = np.insert(self.atom_positions,loc,pos,0)
                self.atom_velocities = np.insert(self.atom_velocities,loc,vel,0)
                self.atom_names.insert(loc,self.atom_names[loc])
                self.atom_types[type] += 1
            else:
                self.atom_positions = np.append(self.atom_positions,pos,0)
                self.atom_velocities = np.append(self.atom_velocities,vel,0)
                self.num_atom_types += 1
                self.atom_names.append(None)
                self.atom_types.append(1)

    def delete_atom(self,index):
        """
        Deletes the atom number index.
        """
        type = 0
        loc = 0
        for i in range(len(self.atom_types)):
            loc += self.atom_types[i]
            type = i
            if loc > index:
                break
        self.num_atoms -= 1
        self.atom_positions = np.delete(self.atom_positions,index,0)
        self.atom_velocities = np.delete(self.atom_velocities,index,0)
        self.atom_names.pop(index)
        if self.atom_types[type] == 1:
            self.atom_types.pop(type)
            self.num_atom_types -= 1
        else:
            self.atom_types[type] -= 1

    def change_atom_type(self,index,type):
        """
        Changes the type of atom number index to the type given.
        """
        pos = self.atom_positions[index]
        vel = self.atom_velocities[index]
        convention = self.convention
        self.delete_atom(index)
        self.add_atom(pos,type,convention,vel)

    def sort_atoms(self):
        """
        Sorts the atoms into alphabetical order by their names. Also updates the
        atomic positions and velocities.
        """
        def argsort(seq):
            return sorted(range(len(seq)),key=seq.__getitem__)
        index = argsort(self.atom_names)
        print index
        print len(index)
        print self.num_atoms
        print self.atom_names
        print self.atom_positions
        self.atom_names.sort()
        temp_atoms = np.zeros_like(self.atom_positions)
        temp_velocities = np.zeros_like(self.atom_positions)
        for i in range(len(index)):
            temp_atoms[i] = self.atom_positions[index[i]]
            temp_velocities[i] = self.atom_velocities[index[i]]
        self.atom_positions = temp_atoms
        self.atom_velocities = temp_velocities

    def lin_trans(self,matrix):
        """
        Applies the linear transformation in numpy array, matrix, to the unit
        cell vectors. Atomic positions are unaffected by the transformation.
        """
        old_conv = self.convention
        self.convention = "Direct"
        new_vec = np.dot(matrix,self.cell_vec.transpose()).transpose()
        self.cell_vec = new_vec
        self.convention = old_conv

    def simple_supercell(self,a,b,c):
        """
        Creates a supercell of the unit cell with lattice vectors scaled by the
        three integers a, b, and c. If a, b, or c, are not integers, they are 
        rounded to the nearest whole number.
        """
        a = int(round(a)); b = int(round(b)); c = int(round(c))
        self.convetion = "Direct"
        list_ = np.array([a,b,c])
        # Create new lattice vectors
        self.cell_vec = (list_*self.cell_vec.transpose()).transpose()
        # Create new atomic positions and update corresponding atom names
        new_positions = np.array(self.atom_positions)
        new_velocities = np.array(self.atom_velocities)
        new_atom_names = self.atom_names
        for i in range(len(list_)):
            temp_positions = np.array(new_positions)
            temp_velocities = np.array(new_velocities)
            temp_names = new_atom_names
            for j in range(1,list_[i]):
                temp_positions[:,i] = temp_positions[:,i] + float(j)
                new_positions = np.append(new_positions,temp_positions,axis=0)
                new_velocities = np.append(new_velocities,temp_velocities,axis=0)
                new_atom_names.extend(temp_names)
            new_positions[:,i] = new_positions[:,i]/float(list_[i])
        self.atom_positions = new_positions
        self.atom_velocities = new_velocities
        self.num_atoms = len(new_positions)
        self._atom_names = new_atom_names
        # Sort the atoms by thier type.
        self.sort_atoms()
        # Update the number of atoms of each type.
        #self.num_atom_types = len(set(self.atom_names))  # Unneccessary
        self.atom_types = []
        for i in range(self.num_atom_types):
            self.atom_types.append(
                    self._atom_names.count(list(set(self._atom_names))[i]))

    def supercell(self,matrix):
        """
        Creates a supercell of the unit cell with lattice vectors scaled by
        matrix. Atoms are then populated in the supercell according to the
        continued lattice. Matrix can either be a 3x3 or 1x3 array-like
        object.
        """
        # Next three lines are probably unneccessary as long as matrix can be
        # turned into a numpy array.
        #if isinstance(matrix,np.ndarray):
        #    matrix = matrix.flatten()
        #else:
        matrix = np.array(matrix).flatten()
        # Check if matrix has only 3 entries, and if so, convert it to the full
        # matrix required in the general case.
        if len(matrix) == 3:
            temp = np.identity(3)
            for i in range(3):
                temp[i,i] = matrix[i]
            matrix = temp.flatten()

    def recip_lat(self):
        """
        Calculates the reciprocal lattice vectors of the unit cell vectors.
        The vectors are returned as rows of a numpy array. The factor of 2*pi 
        is included in the reciprocal vectors.
        """
        self.set_scale()
        recip = 2.0*np.pi*np.linalg.inv(self.cell_vec.transpose())
        return recip

    def center_of_mass(self,potcar="POTCAR",conv="direct",vel=False):
        """
        Calculates the center of mass of the unit cell, given atomic masses from
        a POTCAR formatted file. Center of mass can be returned in either
        cartesian coordinates or in terms of the direct lattice of the cell.
        This is specified by conv - if conv starts with C, c, K, or k, then the
        center of mass will be returned in cartesian coords. and will be
        returned in direct coordinates otherwise.
        """
        # Read in atom masses from POTCAR-formatted file.
        reg_str = r"POMASS *= *([0-9]*[.][0-9]*);"
        regex = re.compile(reg_str)
        potfile = open(potcar,"r")
        lines = potfile.read()
        potfile.close()
        masses = np.array(regex.findall(lines),dtype="float")
        if len(masses) != self.num_atom_types:
            print "The number of atomic masses in the file "+potcar
            print "does not match the number of atom types in the unit cell!"
            print "Program will now exit!"
            sys.exit(1)
        # Calculate center of mass of unit cell.
        self.in_cell()
        self.set_convention(conv)
        n = 0
        r_com = np.zeros(3); tot_mass = 0.0
        v_com = np.zeros(3)
        for i in range(self.num_atom_types):
            for j in range(self.atom_types[i]):
                tot_mass += masses[i]
                r_com += masses[i]*self.atom_positions[n]
                v_com += masses[i]*self.atom_velocities[n]
                n += 1
        r_com = r_com/tot_mass
        v_com = v_com/tot_mass
        if vel:
            return r_com,v_com
        else:
            return r_com

    def spacegroup(self,format_=None):
        """
        Calculate the spacegroup of a crystal using the program findsym. The
        spacegroup can be returned as the International Tables number, the
        Hermann-Mauguin symbol, the Schoenflies symbol, or all three (default).
        If the input argument is not None, then variations on "Number",
        "Schoenflies", and "Hermann-Mauguin" will return the corresponding 
        symbols.
        """
        # Create findsym-formatted output string.
        struc_text = self.output_findsym()
        # Call findsym with above string as input.
        args = "findsym"
        findsym = subprocess.Popen(args,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
        output = findsym.communicate(struc_text)
        # Check that findsym ran correctly.
        if output[1] != None:
            print "findsym exited with error:"
            print output[1]
            sys.exit(1)
        os.remove("findsym.log")  # Log file removed if there is not an error.
        # Grep "Space Group" from output text.
        reg_spg = re.compile(r"Space Group *([0-9]*) *([^ ]*) *([^ ]*)")
        results = reg_spg.search(output[0])
        # (Other output text is still hanging around if you want to use it!)
        # Store space group number, Schoenflies symbol, and Hermann-Mauguin
        # symbol.
        number = results.group(1)
        schoenflies = results.group(2)
        hermann_mauguin = results.group(3)
        full = number+" "+schoenflies+" "+hermann_mauguin
        #full += "\n"+a+" "+b+" "+c
        # Return desired space group representation.
        if ["Number","number","N","n"].count(format_):
            return number
        elif ["Schoenflies","schoenflies","schoen","flies","Sch","sch",
                "S","s"].count(format_):
            return schoenflies
        elif ["Hermann-Mauguin","hermann-mauguin","hermann_mauguin",
                "hermannmauguin","Hermann_Mauguin","HermannMauguin","H-M","h-m",
                "HM","hm","H","h","M","m"].count(format_):
            return hermann_mauguin
        else:
            return full

    def generate_z_list(self,tol=1e-8):
        """
        Returns a 2-d list of each plane of atoms, where the first index of list
        contains each plane of atoms along the z-axis, and the second index
        contains a list of atoms in each of the planes.
        """
        z_positions = self.atom_positions[:,2]
        z_sorted = np.argsort(z_positions)
        result = []
        result.append([z_sorted[0]])
        for i in range(1,len(z_sorted)):
            if abs(z_positions[z_sorted[i]] - z_positions[z_sorted[i-1]]) < tol:
                result[-1].append(z_sorted[i])
            else:
                result.append([z_sorted[i]])
        return np.array(result)

    def read_poscar(self,input_,vel=True,vasp5=False):
        # Currently does not support:
        #   - Magnetic Moments
        #   - DFT+U (does that show up in POSCAR/INCAR for each atom?)
        #   - Selective dynamics
        #   - Vasp 5
        """
        Read in unit cell data from a POSCAR like file line by line.
        """
        # Read in the structure name
        self.name = input_.readline()
        # Read in the scale factor
        self._scale = float(input_.readline().split()[0])
        # Read in the unit cell vectors
        self.cell_vec = np.zeros((3,3))
        for i in range(3):
            line = input_.readline().split()
            for j in range(3):
                self.cell_vec[i,j] = float(line[j])
        # Read in the number of atom types, the number of atoms of each type,
        # and the total number of atoms.
        line = input_.readline().split()
        self.num_atom_types = len(line)
        self.atom_types = []
        self._atom_type_names = []
        self.num_atoms = 0
        # Check to see if file is VASP 5 formatted
        try:
            int(line[0])
        except ValueError:
            # file is formatted for VASP 5
            vasp5 = True
            for i in range(len(line)):
                self._atom_type_names.append(line[i])
            line = input_.readline().split()
        for i in range(len(line)):
            self.atom_types.append(int(line[i]))
            self.num_atoms += int(line[i])
        # Read in atomic coordinate convention
        self._convention = input_.readline().split()[0]
        # Read in the atomic positions.
        self.atom_positions = np.zeros((self.num_atoms,3))
        self._atom_names = []
        for i in range(self.num_atoms):
            line = input_.readline().split()
            for j in range(3):
                self.atom_positions[i,j] = float(line[j])
            # Read in the atomic symbol, if it is in the file.
            if len(line) == 4:
                self._atom_names.append(str(line[3]))
            else:
                self._atom_names.append(None)
        # Check for MD velocities.
        self.atom_velocities = np.zeros((self.num_atoms,3))
        line = input_.readline()
        if vel == True and line:
            # Assume there are velocities if file continues after atoms.
            if ['C','c','K','k','\n','\t',' '].count(line[0]):
                self._vel_convention = 'Cartesian'
            else:
                self._vel_convention = 'Direct'
            for i in range(self.num_atoms):
                line = input_.readline().split()
                try:
                    for j in range(3):
                        self.atom_velocities[i,j] = float(line[j])
                except IndexError:
                    break
        else:
            self._vel_convention = 'Cartesian'
            self.atom_velocities = np.zeros((self.num_atoms,3))
        #input_.close()

    def read_ezvasp(self,input_):
        """
        Read in unit cell information from an ezvasp-formatted file.
        """
        pass

    def read_atat(self,input_):
        """
        Read in unit cell information from an ATAT-formatted file.
        """
        self.name = "ATAT_file"
        self._scale = 1.0
        # Read in coordinate vectors
        coords = np.zeros((3,3))
        for i in range(3):
            line = input_.readline().split()
            # Check if coordinate system is given as 3x3 vectors
            # or a,b,c,alpha,beta,gamma
            if len(line) == 6:
                a = float(line[0]) 
                b = float(line[1]) 
                c = float(line[2])
                alpha = float(line[3]) 
                beta = float(line[4])
                gamma = float(line[5])
                coords = six_to_nine(a,b,c,alpha,beta,gamma)
                break
            for j in range(3):
                coords[i,j] = float(line[j])
        # Read in unit cell vectors
        unit_vec = np.zeros((3,3))
        for i in range(3):
            line = input_.readline().split()
            for j in range(3):
                unit_vec[i,j] = float(line[j])
        self.cell_vec = np.dot(coords,unit_vec)
        # Read in the atomic positions
        lines = input_.readlines()
        atom_positions = []
        atom_names = []
        self._convention = 'Cartesian'
        for i,line in enumerate(lines):
            line = line.split()
            if len(line) == 4:
                atom_names.append(str(line.pop(-1)))
                atom_pos = np.array( [ float(i) for i in line ] )
                atom_pos = np.dot(coords,atom_pos)
                atom_positions.append(atom_pos)
            else:
                break
        # Sort atomic positions by atom name
        atom_names = np.array(atom_names)
        atom_positions = np.array(atom_positions)
        sort_pos = np.argsort(atom_names)
        self._atom_names = atom_names[sort_pos]
        self.atom_positions = atom_positions[sort_pos]
        self.num_atoms = len(self.atom_names)
        # Determine distinct atom types
        self.atom_type_names = []
        self.atom_types = []
        self.atom_type_names.append(self.atom_names[0])
        self.atom_types.append(1)
        self.num_atom_types = 1
        for name in self.atom_names[1:]:
            if name == self.atom_type_names[-1]:
                self.atom_types[-1] += 1
            else:
                self.num_atom_types += 1
                self.atom_type_names.append(name)
                self.atom_types.append(1)
        # Set default values for undefined UnitCell parameters
        self._vel_convention = 'Cartesian'
        self.atom_velocities = np.zeros_like(self.atom_positions)

    def read_gulp_input(self,input_):
        """
        Read in unit cell information from a gulp input-formatted file.
        """
        pass
        #input_.readline()  # Discard line containing calculations to run.
        #line = input_.readline()
        #if line.split
        lines = input_.read()
        for i,line in enumerate(lines):
            if line.startswith("cell"):
                a = float(lines[i+1].split()[0])
                b = float(lines[i+1].split()[1])
                c = float(lines[i+1].split()[2])
                alpha = float(lines[i+1].split()[0])
                beta = float(lines[i+1].split()[0])
                gamma = float(lines[i+1].split()[0])
                self.cell_vec = six_to_nine(a,b,c,alpha,beta,gamma)
            elif line.startswith("frac"):
                pass

    def read_gulp_output(self,input_,top=False):
        """
        Read gulp output file to grab atomic coordinates. If top = False, the
        final relaxed positions are taken. If top = True, or no final positions
        exist, then the initial atomic positions in the gulp output file will be
        read.
        """
        # Read in gulp output file
        lines = input_.read()
        # Check if final coordinates are in gulp output file.
        reg_opti = re.compile(r"optimise *- *perform optimisation run")
        if not reg_opti.search(lines):
            top = True
        # Compile regular expressions for getting unit cell data.
        reg_num_atoms = re.compile(r"Total number atoms/shells *= *([0-9]*)")
        if not top:
            reg_a = re.compile(r"\ba *([0-9]*[.][0-9]*) *Angstrom")
            reg_b = re.compile(r"\bb *([0-9]*[.][0-9]*) *Angstrom")
            reg_c = re.compile(r"\bc *([0-9]*[.][0-9]*) *Angstrom")
            reg_alpha = re.compile(r"alpha *([0-9]*[.][0-9]*) *Degrees")
            reg_beta = re.compile(r"beta *([0-9]*[.][0-9]*) *Degrees")
            reg_gamma = re.compile(r"gamma *([0-9]*[.][0-9]*) *Degrees")
            reg_start = re.compile(r"Label *\(Frac\) *\(Frac\) *\(Frac\) *\(Angs\) *\n-*")
            reg_stop = re.compile(r"-*\n\n *Final Cartesian lattice vectors")
            self.name = "Gulp Output Final Relax\n"
        else:
            reg_a = re.compile(r"\ba *= *([0-9]*[.][0-9]*)")
            reg_b = re.compile(r"\bb *= *([0-9]*[.][0-9]*)")
            reg_c = re.compile(r"\bc *= *([0-9]*[.][0-9]*)")
            reg_alpha = re.compile(r"alpha *= *([0-9]*[.][0-9]*)")
            reg_beta = re.compile(r"beta *= *([0-9]*[.][0-9]*)")
            reg_gamma = re.compile(r"gamma *= *([0-9]*[.][0-9]*)")
            reg_start = re.compile(
               r"Label *\(Frac\) *\(Frac\) *\(Frac\) *\(e\) *\(Frac\) *\n-*\n")
            reg_stop = re.compile(r"-*\n*\**\n\* *General input information")
            self.name = "Gulp Output Initial Geometry\n"
        # Input unit cell vectors.
        a  = float(reg_a.search(lines).group(1))
        b  = float(reg_b.search(lines).group(1))
        c  = float(reg_c.search(lines).group(1))
        alpha = float(reg_alpha.search(lines).group(1))
        beta = float(reg_beta.search(lines).group(1))
        gamma = float(reg_gamma.search(lines).group(1))
        self.cell_vec = six_to_nine(a,b,c,alpha,beta,gamma)
        # Input other unit cell attributes.
        self._scale = 1.0
        self._convention = "Direct"
        self.num_atoms = int(reg_num_atoms.search(lines).group(1))
        # Input atomic positions.
        start = reg_start.search(lines).end()
        stop = reg_stop.search(lines).start()
        atom_lines = lines[start+1:stop-1].split('\n')
        self._atom_names = []
        self.atom_positions = np.zeros((self.num_atoms,3))
        for i in range(len(atom_lines)):
            line = atom_lines[i].split()
            self._atom_names.append(line[1])
            if len(line) <= 8:
                self.atom_positions[i,0] = float(line[3])
                self.atom_positions[i,1] = float(line[4])
                self.atom_positions[i,2] = float(line[5])
            else:
                self.atom_positions[i,0] = float(line[3])
                self.atom_positions[i,1] = float(line[5])
                self.atom_positions[i,2] = float(line[7])
        # Get atom type information from atomic positions.
        self.num_atom_types = len(set(self._atom_names))
        self.atom_types = []
        for i in range(self.num_atom_types):
            self.atom_types.append(
                    self._atom_names.count(list(set(self._atom_names))[i]))
        self._vel_convention = 'Direct'
        self.atom_velocities = np.zeros((self.num_atoms,3))

    def read_lammps_dump(self,input_):
        # Currently does not support:
        #   - non-orthogonal cells
        #   - MD velocities
        """
        Read LAMMPS dump file in to grab the last set of atomic coordinates.
        """
        self.name = "LAMMPS Dump file.\n"
        self._scale = float(1.0)
        # Currently, assume that input_ is a string containing only the last
        # set of dump data.
        lines = input_.split('\n')
        for i in range(len(lines)):
            if lines[i] == "ITEM: NUMBER OF ATOMS":
                self.num_atoms = int(lines[i+1])
                break
        # Currently only works for orthogonal unit cells!
        self.cell_vec = np.zeros((3,3))
        for i in range(len(lines)):
            if lines[i] == "ITEM: BOX BOUNDS":
                self.cell_vec[0,0] = lines[i+1].split()[1]
                self.cell_vec[1,1] = lines[i+2].split()[1]
                self.cell_vec[2,2] = lines[i+3].split()[1]
                break
        self.atom_positions = np.zeros((self.num_atoms,3))
        self._atom_names = []
        self._convention = "Direct"
        flag = False; j = 0
        for i in range(len(lines)):
            if flag == True:
                self._atom_names.append(lines[i].split()[1])
                self.atom_positions[j,0] = float(lines[i].split()[2])
                self.atom_positions[j,1] = float(lines[i].split()[3])
                self.atom_positions[j,2] = float(lines[i].split()[4])
                j += 1
            if lines[i] == "ITEM: ATOMS id type xs ys zs":
                flag = True
                j = 0
            if j == self.num_atoms:
                break
        self.num_atom_types = len(set(self._atom_names))
        self.atom_types = []
        for i in range(self.num_atom_types):
            self.atom_types.append(
                    self._atom_names.count(list(set(self._atom_names))[i]))
        # Currently does not read in atomic velocities
        self._vel_convention = 'Direct'
        self.atom_velocities = np.zeros((self.num_atoms,3))

    def output_vasp(self):
        # Currently does not support:
        #   - Magnetic moments
        #   - DFT+U
        #   - Vasp 5
        """
        Returns a string containing the unit cell information in a
        VASP-formatted file.
        """
        if self.name[-1] == "\n":
            output = self.name
        else:
            output = self.name+"\n"
        output += str(self.scale)+"\n"
        # Output unit cell vectors
        for i in range(3):
            for j in range(3):
                output += str(self.cell_vec[i,j])+" "
            output += "\n"
        # Output numbers and types of atoms
        for i in range(len(self.atom_types)):
            output += str(self.atom_types[i])+" "
        output += "\n"
        output += self.convention+"\n"
        # Output atomic positions
        for i in range(self.num_atoms):
            for j in range(3):
                output += str(self.atom_positions[i,j])+" "
            if self.atom_names[i] != None:
                output += self.atom_names[i]
            output += "\n"
        # Output atomic velocities, if they are non-zero
        if self.atom_velocities.any():
            output += self.vel_convention+"\n"
            for i in range(self.num_atoms):
                for j in range(3):
                    output += str(self.atom_velocities[i,j])+" "
                output += "\n"
        return output

    def output_ezvasp(self):
        # Currently does not support:
        #   - MD velocities
        #   - Vasp 5?
        #   - Magnetic moments
        #   - DFT+U?
        """
        Returns a string containing the unit cell information in ezvasp
        formatting, ready to be appended to the end of an existing vasp.in-like
        file.
        """
        # Create ezvasp-formatted text
        if self.name[-1] == "\n":
            output = self.name
        else:
            output = self.name+"\n"
        output += str(self.scale)+"\n"
        # Output unit cell vectors
        for i in range(3):
            for j in range(3):
                output += str(self.cell_vec[i,j])+" "
            output += "\n"
        output += self.convention+"\n"
        # Output atomic positions
        for i in range(self.num_atoms):
            for j in range(3):
                output += str(self.atom_positions[i,j])+" "
            if self.atom_names[i] != None:
                output += self.atom_names[i]
            output += "\n"
        return output

    def output_gulp(self):
        """
        Returns a string containing the unit cell information in a
        GULP-formatted file.
        """
        self.convention = 'Direct'
        self.set_scale()
        calculations = "opti conp prop"
        output = calculations+"\n"
        if self.name[-1] == "\n":
            output += "title\n"+self.name+"end\n"
        else:
            output += "title\n"+self.name+"\nend\n"
        # Output unit cell vectors
        output += "vector\n"
        for i in range(3):
            for j in range(3):
                output += str(self.cell_vec[i,j])+" "
            output += "\n"
        # Output atomic positions in direct coordinates
        output += "frac\n"
        for i in range(self.num_atoms):
            output += self.atom_names[i]+" "
            for j in range(3):
                output += str(self.atom_positions[i,j])+" "
            output += "\n"
        # Default and non-useful interatomic potential data
        output += "species\n"
        output += "library\n"
        output += "output movie xyz cluster\n"
        return output

    def output_lammps(self):
        # Currently does not support:
        #   - MD velocities
        """
        Returns a string containing the unit cell information in a
        LAMMPS-formatted file. Only works for orthogonal unit cells!
        """
        self.set_scale()
        self.to_cart()
        if self.name[-1] == "\n":
            output = self.name+"\n"
        else:
            output = self.name+"\n\n"
        output += str(self.num_atoms)+" atoms\n\n"
        output += str(self.num_atom_types)+" atom types\n\n"
        #unit cell geometry goes here. Figure out a better way later.
        a = np.linalg.norm(self.cell_vec[0])
        b = np.linalg.norm(self.cell_vec[1])
        c = np.linalg.norm(self.cell_vec[2])
        output += "0 "+str(a)+" xlo xhi\n"
        output += "0 "+str(b)+" ylo yhi\n"
        output += "0 "+str(c)+" zlo zhi\n\n"
        # Output atomic types and positions
        output += "Atoms\n\n"
        j=0; n=1
        for i in range(self.num_atoms):
            if j >= self.atom_types[n-1]:
                n += 1
                j = 0
            output += str(i+1)+" "+str(n)
            for k in range(3):
                output += " "+str(self.atom_positions[i,k])
            output += "\n"
            j += 1
        return output

    def output_findsym(self):
        """
        Returns a string containing the unit cell information in a
        findsym-formatted file.
        """
        self.set_scale()
        self.convention = 'Direct'
        if self.name[-1] == "\n":
            output = self.name
        else:
            output = self.name+"\n"
        output += "0\n"  # Use default symmetry tolerance (1e-6)
        output += "1\n"  # Set lattice parameter format to vectors
        # Output unit cell vectors as "lattice parameter vectors"
        for i in range(3):
            for j in range(3):
                output += str(self.cell_vec[i,j])+" "
            output += "\n"
        output += "1\n"  # Set unit cell vector format to vectors
        # Unit cell vectors chosen as identity matrix - all relevant information
        # is contained in the "lattice parameter vectors"
        output += "1.0 0.0 0.0\n"
        output += "0.0 1.0 0.0\n"
        output += "0.0 0.0 1.0\n"
        output += str(self.num_atoms)+"\n"
        # Output list of atomic type integers
        n = 1
        for i in range(self.num_atom_types):
            for j in range(self.atom_types[i]):
                output += str(n)+" "
            n += 1
        output += "\n"
        # Output atomic positions
        for i in range(self.num_atoms):
            for j in range(3):
                output += str(self.atom_positions[i,j])+" "
            output += "\n"
        return output

def cart_to_direct(lattice,property):
    """
    Converts the 1- or 2-D array property from cartesian to direct coordinates
    using the 3x3 array lattice. The inputs lattice and property are both
    assumed to be numpy arrays. Property is assumed to be an Nx3 array
    containing N vectors in cartesian coordinates.
    """
    return np.dot(property,np.linalg.inv(lattice))

def direct_to_cart(lattice,property):
    """
    Converts the 1- or 2-D array property from direct to cartesian coordinates
    using the 3x3 array lattice. The inputs lattice and property are both
    assumed to be numpy arrays. Property is assumed to be an Nx3 array
    containing N vectors in direct coordinates.
    """
    return np.dot(property,lattice)

def six_to_nine(a,b,c,alpha,beta,gamma):
    """
    Converts unit cell parameters a,b,c,alpha,beta,gamma to lattice vectors.
    Angles should be given in degrees and distances in angstroms.
    """
    alpha = alpha*2*np.pi/360.
    beta = beta*2*np.pi/360.
    gamma = gamma*2*np.pi/360.
    lat = np.zeros((3,3))
    lat[0,:] = np.array([a,0,0])
    lat[1,:] = np.array([b*np.cos(gamma),b*np.sin(gamma),0])
    lat[2,0] = c*np.cos(beta)
    lat[2,1] = c*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)
    lat[2,2] = c*np.sqrt(1-(np.cos(alpha)**2 + np.cos(beta)**2 
                - 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))/np.sin(gamma)**2)
    lat = np.round(lat,decimals=8)
    return lat

def displacements(cell_1,cell_2,conv="C",flag=False):
    """
    Calculate the displacements of atoms in unit cell 2 from the atomic
    coordinates of unit cell 1. If flag is False (default), return the
    displacement vectors of each atom. If flag is True, return the
    magnitude of the displacement and the [uvw] direction of displacement.
    """

    def gcd(a,b):
        """
        Calculate the greatest common divisor of 2 numbers.
        """
        while b !=0:
            (a,b) = (b,a%b)
        return a

    def gcdd(*args):
        """
        Calculate the greatest common divisor of a set of numbers.
        """
        return reduce(gcd,args)

    # Make sure that the unit cell vectors are identical in the two cells.
    # and that they have the same number of atoms.
    #   I'll put that in later!
    if cell_1.num_atoms == cell_2.num_atoms:
        atoms = cell_1.num_atoms
    else:
        print "Unit cells have different numbers of atoms!"
        return
    # Loop over each atom and calculate the displacement.
    cees = ['c','C','k','K']
    cell_1.convention = 'Direct'
    cell_2.convention = 'Direct'
    cell_1.scale = 1.0
    cell_2.scale = 1.0
    disp = np.zeros((atoms,3))
    if flag == True:
        mag = np.zeros(atoms)
        uvw = np.zeros((atoms,3))
    for i in range(atoms):
        disp[i] = cell_2.atom_positions[i] - cell_1.atom_positions[i]
        for j in range(3):
            if disp[i,j] > 0.5:
                disp[i,j] = (cell_2.atom_positions[i,j] - 1 -
                cell_1.atom_positions[i,j])
        if cees.count(conv[0]):
            disp[i] = np.dot(cell_1.cell_vec.transpose(),disp[i])
        if flag == True:
            mag[i] = np.linalg.norm(disp[i])
            uvw[i] = disp[i]/abs(gcdd(disp[i,0],disp[i,1],disp[i,2]))
    if flag == False:
        return disp
    else:
        return mag,uvw

if __name__ == "__main__":
    poscar = open(sys.argv[1],'r')
    unit_cell = UnitCell(poscar)
    outfile = unit_cell.spacegroup()
    print outfile
    sys.exit()
