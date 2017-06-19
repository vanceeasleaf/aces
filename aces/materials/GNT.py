# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2017-06-19 13:16:44
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-19 13:17:00
from aces.materials import Material
from math import pi, cos, sin
from aces.materials.graphene import structure as graphene
from aces.materials import atomic
import numpy as np


class structure(Material):

    def set_parameters(self):
        self.type = 'zigzag'
        self.phi = 2 * pi

    def setup(self):
        self.enforceThick = False

    def lmp_structure(self):
        atoms = graphene(
            dict(
                latx=self.latx,
                laty=self.laty,
                latz=self.latz,
                gnrtype=self.type)).lmp_structure()
        atomic.center_box(atoms)
        self.radius = atoms.cell[1][1] / self.phi
        for atom in atoms:
            atom.position = self.trans(atom.position, r=self.radius)

        atoms.center(vacuum=10, axis=[1, 2])
        atoms.center(axis=[0])
        return atoms

    def trans(self, pos, r=1):
        x, y, z = np.array(pos)
        x1 = x
        t = y / r
        y1 = sin(t) * (r + z)
        z1 = r - cos(t) * (r + z)
        return np.array([x1, y1, z1])
