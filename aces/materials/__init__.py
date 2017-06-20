# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2017-06-17 13:34:52
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-20 16:08:25


import aces.tools as tl
from ase import Atoms
from aces import default
from aces.io.vasp import writevasp
from ase import io
from aces.tools.Units import Units
from aces import config
from ase.dft.kpoints import ibz_points
from aces.io.lammps.lammpsdata import lammpsdata
import numpy as np
from aces.env import PROJHOME
import atomic


class Material:

    def __init__(self, opt={}):

        # all the values needed
        # unit might be changed by opt but need to be used first
        self.__dict__ = dict(self.__dict__,
                             **default.default)
        if 'units' in opt:
            self.units = opt['units']
        self.units = Units(self.units)
        self.elements = ['C', 'N', 'B']
        self.set_parameters()
        self.__dict__ = dict(self.__dict__, **opt)
        self.super_setup()

    def __getattr__(self, attr):
        if attr == "dim":
            return tl.toString(self.supercell)
        if attr == "cores":
            return int(self.nodes) * int(self.procs)

    # to be overided
    def set_parameters(self):
        pass

    def super_setup(self):
        self.units = Units(self.units)
        self.prepare_lammps()

        self.prepare_phonts()
        self.bandpoints = ibz_points['fcc']
        self.bandpoints['X'] = [.5, 0, 0]
        self.bandpoints['Y'] = [0, 0.5, 0]
        self.bandpath = ['Gamma', 'X', 'Y', 'Gamma']
        self.premitive = np.eye(3)
        if not self.useS3:
            self.supercell3 = self.supercell
        self.setup()
        if self.atomfile:
            atoms = io.read(
                str(PROJHOME + "/data/" + self.atomfile), format="vasp")
            self.atoms = atoms.repeat([self.latx, self.laty, self.latz])
            self.atoms.center()
        else:
            self.atoms = self.lmp_structure()
        if self.dimension == 1:
            self.masses += "\nfix   1d all setforce NULL 0. 0.\nvelocity  all set NULL 0.0 0.0 units box"
        elif self.dimension == 2:
            self.masses += "\nfix   1d all setforce NULL NULL 0.\nvelocity  all set NULL NULL 0.0 units box"
    # to be overided

    def setup(self):
        pass

    def prepare_lammps(self):
        self.potential = 'pair_style	tersoff\npair_coeff	* * %s/BNC.tersoff  %s' % (
            config.lammpspot, ' '.join(self.elements))
        self.dump = "dump_modify dump1 element %s" % (' '.join(self.elements))
        masses = atomic.getMassFromLabel(self.elements)
        self.masses = '\n'.join(["mass %d %f" % (i + 1, mass)
                                 for i, mass in enumerate(masses)])
        m = self
        units = self.units
        m.kb = units.boltz
        m.nktv = units.nktv2p
        if(m.method == "nvt"):
            m.xp = 0
        m.dtime = m.timestep * 100
        m.tcfactor = units.tcfactor
        m.excNum = m.aveRate / m.excRate
        m.swapEnergyRate = m.swapEnergy / (m.excRate * m.timestep)

    def prepare_phonts(self):
        masses = atomic.getMassFromLabel(self.elements)
        self.phontsmasses = '\n'.join(
            ["%s %f 0.0" % (label, mass) for label, mass in zip(self.elements, masses)])

    def structure(self):

        self.write()

    # to be overrided
    def lmp_structure(self):
        atoms = Atoms()
        return atoms

    def write(self):
        self.watoms(self.atoms)

    def watoms(self, atoms):
        atoms.write("structure.xyz")
        writevasp(atoms)
        # write_vasp("POSCAR",atoms,sort="True",direct=True,vasp5=True)
        self.POSCAR2data()
        if len(atoms) < 1000:
            atoms.write('structure.png')

    def writeatoms(self, atoms, label='atoms'):
        tl.mkcd(label)
        self.watoms(atoms)
        tl.cd('..')

    def getatomicstyle(self):
        a = "atom_style atomic"
        if self.creatbonds > 0.0:
            a = "atom_style bond\natom_modify sort 0 1.\ncomm_modify  cutoff 2.0 "
        return a

    def POSCAR2data(self):
        atoms = io.read('POSCAR')
        m = self
        atoms.set_pbc([m.xp, m.yp, m.zp])
        # debug(atoms.cell)
        a = lammpsdata(atoms, self.elements)
        rot = a.writedata(filename="structure", creatbonds=self.creatbonds)
        d, p, d1, p1 = rot
        np.savetxt('POSCARrot', np.r_[d, p, d1, p1])
        # debug(rot)
        return rot

    def atoms_from_dump(self, filename):
        from atomic import atoms_from_dump as afd
        atoms = afd(filename=filename, elements=self.elements)
        m = self
        atoms.set_pbc([m.xp, m.yp, m.zp])
        return atoms

    def dump2POSCAR(self, dumpname, poscar='POSCAR', rotate=True):
        atoms = self.atoms_from_dump(dumpname)
        if rotate:
            rot = np.loadtxt(tl.dirname(dumpname) + '/POSCARrot')
            d, p, d1, p1 = rot[:3], rot[3], rot[4:7], rot[7]
            atoms.rotate(d1, -p1, rotate_cell=True)
            atoms.rotate(d, -p, rotate_cell=True)
        # write_vasp(poscar,atoms,sort="True",direct=True,vasp5=True)
        writevasp(atoms, poscar)

    def getboxrange(self):
        file = open("range")
        for i in range(5):
            file.next()
        xlo, xhi = map(float, file.next().split()[:2])
        ylo, yhi = map(float, file.next().split()[:2])
        zlo, zhi = map(float, file.next().split()[:2])
        return (xlo, xhi, ylo, yhi, zlo, zhi)

    def getxrange(self):
        file = open('minimize.xyz')
        n = file.next().split()[0]
        n = int(n)
        file.next()
        xmin = 100000
        xmax = -100000
        ymin = 100000
        ymax = -100000
        zmin = 100000
        zmax = -100000
        for i in range(n):
            label, x, y, z = file.next().split()
            x, y, z = map(float, [x, y, z])
            xmin = min(x, xmin)
            xmax = max(x, xmax)
            ymin = min(y, ymin)
            ymax = max(y, ymax)
            zmin = min(z, zmin)
            zmax = max(z, zmax)
        return (xmin, xmax, ymin, ymax, zmin, zmax)

    def postMini(self):
        xlo, xhi, ylo, yhi, zlo, zhi = self.getboxrange()
        xlo0, xhi0, ylo0, yhi0, zlo0, zhi0 = self.getxrange()
        if(self.xp == 0):
            xlo = xlo0
            xhi = xhi0
        if(self.yp == 0):
            ylo = ylo0
            yhi = yhi0
        if(self.zp == 0):
            zlo = zlo0
            zhi = zhi0
        lx = xhi - xlo
        ly = yhi - ylo
        lz = zhi - zlo
        if(self.enforceThick):
            self.zfactor = lz / self.thick
        else:
            self.zfactor = 1
        self.S = ly * lz
        self.box = (xlo, xhi, ylo, yhi, zlo, zhi, lx, ly, lz)
