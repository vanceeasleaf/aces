from ..lineManager import lineManager
from aces.scanf import sscanf
import numpy as np


class phononyaml:

    def __init__(self, filename):
        lm = self.lm = lineManager(filename)
        self.mesh = sscanf(lm.getLine(0), 'mesh: [    %d,    %d,     %d ]')
        self.nqpoint, = sscanf(lm.getLine(1), 'nqpoint: %d')
        self.natom, = sscanf(lm.getLine(2), 'natom:   %d')
        if 'reci' in lm.getLine(3):
            self.off = 4
        else:
            self.off = 0
        self.nbranch = 3 * self.natom
        self.lenatom = 4
        self.leneigvec = self.natom * self.lenatom + 1
        self.lenbranch = self.leneigvec + 2
        self.lenband = self.lenbranch * self.nbranch + 1
        self.lenqpoint = self.lenband + 3

    def qposition(self, iqp, i=False):
        lm = self.lm
        iline = 4 + self.off + iqp * self.lenqpoint
        if i:
            return iline
        return sscanf(
            lm.getLine(iline), '- q-position: [    %f,    %f,    %f ]')

    def frequency(self, iqp, ibr, i=False):
        lm = self.lm
        iline = self.qposition(iqp, i=True) + 3 + ibr * self.lenbranch + 1
        if i:
            return iline
        freq, = sscanf(lm.getLine(iline), '    frequency:    %f')
        return freq

    def atom(self, iqp, ibr, ia):
        lm = self.lm
        pos = np.zeros(3, dtype=np.complex)
        iline = self.frequency(iqp, ibr, True) + 2 + self.lenatom * ia + 1
        lm.moveto(iline)
        for i in range(3):
            x1, x2 = sscanf(lm.nextLine(), '      - [  %f,  %f ]')
            pos[i] = x1 + 1j * x2
        return pos

    def atoms(self, iqp, ibr):
        pos = np.zeros([self.natom, 3], dtype=np.complex)
        for i in range(self.natom):
            pos[i] = self.atom(iqp, ibr, i)
        return pos
