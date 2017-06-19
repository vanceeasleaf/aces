# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2017-06-19 13:09:14
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-19 13:09:20
from aces.tools import parseyaml
import numpy as np


class phononyaml:

    def __init__(self, filename):
        self.data = parseyaml(filename)
        data = self.data
        self.nqpoint = int(data['nqpoint'])
        self.natom = int(data['natom'])

        self.reciprocal_lattice = np.array(data['reciprocal_lattice'])
        self.phonon = data['phonon']
        self.nbranch = 3 * self.natom

    def qposition(self, iqp):
        return map(float, self.phonon[iqp]['q-position'])

    def band(self, iqp, ibr):
        return self.phonon[iqp]['band'][ibr]

    def frequency(self, iqp, ibr):
        return float(self.band(iqp, ibr)['frequency'])

    def atoms(self, iqp, ibr):
        vec = np.array(self.band(iqp, ibr)['eigenvector'])
        return vec[:, :, 0] + vec[:, :, 1] * 1j
