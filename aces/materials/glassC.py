# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2017-06-12 18:37:34
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-26 14:24:56
#
"""pressed glass carbon synthised in

[1] M. Hu, et al.,
“Compressed glassy carbon: An ultrastrong
 and elastic interpenetrating graphene network,”
 no. June, pp. 1–8, 2017.

and discripted in http://mp.weixin.qq.com/s?
__biz=MzAwNzU5NjY5MA==&mid=2658343029&idx=1&sn=
b5f48fbb3a8ced653770bcbd6587b21b&chksm=80fc0788b78b8e9e74c1
e112946532d0f58ba14cffa8848e64144564d7d1d1f05fe57609d882&
mpshare=1&scene=1&srcid=0612FlqVOY1vYShoZDxFcd6H#rd

the structure can be build as follows:
1. build a ziggag graphene block with 2x1x1 unit
2. find the rightest atoms as rotation center,
 copy and rotate the block along y axis with -120
"""
from aces.materials import Material
from ase import Atoms
import atomic
from aces.materials.graphene import structure as graphene
import numpy as np


class structure(Material):

    def setup(self):
        self.enforceThick = False

    def lmp_structure(self):
        atoms = Atoms()
        unit = graphene(
            dict(
                latx=2,
                laty=1,
                latz=1,
                gnrtype='zigzag')).lmp_structure()
        # del unit[unit.positions[:,0].argsort()[-2:]]
        ly = unit.cell[1][1]
        lidx = unit.positions[:, 0].argsort()[0]
        left_of_unit = unit[lidx].position
        # the center must be assign because the atoms are centered
        unit.rotate('y', 30.0 / 180 * np.pi, center=left_of_unit)
        ridx = unit.positions[:, 0].argsort()[-1]
        right_of_unit = unit[ridx].position
        unit1 = unit.copy()
        unit1.rotate('y', 120.0 / 180 * np.pi, center=right_of_unit)
        ridx1 = unit1.positions[:, 0].argsort()[-1]
        right_of_unit1 = unit1[ridx1].position
        # cell[0,0]
        lx = right_of_unit1[0] - left_of_unit[0]

        unit2 = unit.copy()
        # translate unit2 but don't traslate along y
        deta = right_of_unit1 - right_of_unit + [0, 0, 1.42 * np.sqrt(3) / 2]
        deta[1] = 0
        unit2.translate(deta)
        lidx2 = unit2.positions[:, 0].argsort()[0]
        left_of_unit2 = unit2.positions[lidx2]
        unit3 = unit2.copy()
        unit3.rotate('y', 120.0 / 180 * np.pi, center=left_of_unit2)

        lz = left_of_unit2[2] - right_of_unit[2] + 1.42 * np.sqrt(3) / 2

        atoms.extend(unit)
        atoms.extend(unit1)
        atoms.extend(unit2)
        atoms.extend(unit3)
        atoms.set_cell([lx, ly, lz])
        # important for distance calculation,default all False
        atoms.set_pbc([True] * 3)
        atoms.center()
        atoms = atomic.get_unique_atoms(atoms)

        # prevent xs~=1.0
        atoms.translate([0.01, 0, 0])
        atomic.wrap(atoms)

        return atoms
