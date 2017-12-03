# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2017-06-13 21:41:19
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-10-28 12:03:21
from aces.materials import Material
import aces.tools as tl
from aces import config
from ase import io


class structure(Material):

    def set_parameters(self):
        self.type = 'zigzag'
        self.poscar = '2l1'
        self.elements = ['Si']

    def setup(self):
        self.bandpoints = {
            'K': [0.5, 0.5, 0],
            'Gamma': [0, 0, 0],
            'X': [0.5, 0.0, 0],
            'Y': [0.0, 0.5, 0]
        }
        self.bandpath = ['X', 'Gamma', 'Y']
        tl.debug("graphene potential chosen:Si.tersoff.mod")
        self.potential = 'pair_style	tersoff/mod\n' +\
            'pair_coeff	* * %s/Si.tersoff.mod Si' % (
                config.lammpspot)

    def lmp_structure(self):
        import os
        file = os.path.join(
            os.path.dirname(__file__), 'data', self.poscar + ".POSCAR")
        # because file is unicode in default and ase check if file is not str
        # then look it as a file
        file = str(file)
        assert tl.exists(file)
        print(file)
        atoms = io.read(file)
        if not self.type == "zigzag":
            newatoms = atoms.copy()
            cell = [atoms.cell[1, 1], atoms.cell[0, 0], atoms.cell[2, 2]]
            positions = atoms.positions[:, [1, 0, 2]]
            newatoms.set_cell(cell)
            newatoms.set_positions(positions)
            atoms = newatoms
        atoms = atoms.repeat([self.latx, self.laty, self.latz])
        return atoms
