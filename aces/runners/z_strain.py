# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2017-06-26 14:13:40
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-26 17:25:41

from aces.runners import Runner
from aces.runners.minimize import minimize as minimize_input
import aces.tools as tl
from aces.io.vasp import writevasp
import numpy as np
from ase import io


class runner(Runner):

    def generate(self):
        self.run()

    def creatmini(self, dir):
        cur = tl.pwd()
        tl.mkdir(dir)
        tl.cd(dir)
        minimize_input(self.m)
        tl.cd(cur)

    def run_next(self, dir0, dir, lz):
        m = self.m

        atoms = io.read("POSCAR_" + dir0)
        m.forceThick = True
        cell = atoms.get_cell()
        cell[2][2] = lz
        atoms.set_cell(cell, scale_atoms=True)
        m.atoms = atoms
        self.creatmini(dir)
        atoms = m.atoms_from_dump('%s/range' % dir)
        writevasp(atoms, "POSCAR_" + dir)

    def run(self):
        m = self.m
        dir0 = 'minimize'
        atoms = m.dump2POSCAR('%s/range' % dir0)
        tl.cp("POSCAR", "POSCAR_minimize")
        cell = atoms.get_cell()
        r = np.arange(1.0, 0.5, -0.01) * cell[2][2]
        for i, lz in enumerate(r):
            dir = 'minimize%d' % i
            self.run_next(dir0, dir, lz)
            dir0 = dir
