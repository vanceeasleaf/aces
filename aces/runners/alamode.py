# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2015-11-01 03:24:48
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-19 12:48:32

import aces.tools as tl
from ase.io.vasp import write_vasp
from ase import io
from aces import config
from aces.runners.phonopy import runner as Runner
from aces.f import read_forces, matrixFormat


class runner(Runner):

    def generate(self):
        # self.minimizePOSCAR()
        tl.shell_exec('cp minimize/POSCAR .')
        self.get_almin()
        self.displacements()
        files, self.NDATA = self.getfiles()
        self.getvasprun(files)
        self.getdispforce(files)
        self.get_fitin()
        self.getfcs()
        self.get_anphonin()
        self.run_anphonin()

    def run_anphonin(self):
        tl.passthru(config.mpirun + str(self.m.cores) + config.anphon +
                    " band.in > band.out")
        tl.passthru(config.mpirun + str(self.m.cores) + config.anphon +
                    " dos.in > dos.out")
        tl.passthru(config.mpirun + str(self.m.cores) + config.anphon +
                    " tc.in > tc.out")

    def getfcs(self):
        tl.passthru(config.alm + "< fit.in > fit.out")
        assert tl.exists("alm.fcs")
        assert tl.exists("alm.xml")

    def getfiles(self):
        files = tl.shell_exec("ls *.POSCAR|sort -n").split('\n')
        assert len(files) > 0 and not files[0] == ''
        return files, len(files)

    def get_fitin(self):
        content = tl.read('alm.in').replace('suggest', 'fitting')
        fitting = """&fitting
        \tNDATA = %d
        \tDFILE = disp_all.dat
        \tFFILE = force_all.dat
        /
        """ % self.NDATA
        fitting = tl.trimhead(fitting)
        tl.write(content + fitting, 'fit.in')

    def getdispforce(self, files):
        force = ""
        disp = ""
        orig = io.read('POSCAR-supercell').positions
        for dir0 in files:
            forcearr = read_forces(
                'dirs/dir_%s/vasprun.xml' % dir0) / 25.7110  # in Rd/bohr
            force += matrixFormat(forcearr)
            disparr = (io.read('dirs/dir_%s/POSCAR' % dir0).positions - orig
                       ) * 1.889726  # in bohr
            disp += matrixFormat(disparr)
        tl.write(disp, "disp_all.dat")
        tl.write(force, "force_all.dat")

    def displacements(self):
        tl.passthru(config.alm + "< alm.in > alm.out")
        files = tl.shell_exec("ls *pattern*").split()
        tl.passthru(config.almdisp + self.m.toString(files))

    def get_almin(self):
        m = self.m
        m.atoms = io.read('POSCAR')
        atoms = m.atoms.repeat(m.supercell)

        general = """&general
        \tPREFIX = alm
        \tMODE = suggest
        \tNAT = %s; NKD = %s
        \tKD = %s
        /
        """ % (len(atoms), len(m.elements), m.toString(m.elements))

        interaction = """&interaction
        \tNORDER = 2
        /
        """
        cell = """&cell
        \t1.889726
        \t%s
        /
        """ % ('\n  '.join([m.toString(atoms.cell[i]) for i in range(3)]))

        cutoff = """&cutoff
        \t*-* None %f
        /
        """ % (self.m.shengcut * 1.889726)

        pos = '  \n'.join([
            '\t%s ' % (m.elements.index(a.symbol) + 1) +
            m.toString(atoms.get_scaled_positions()[i])
            for i, a in enumerate(atoms)
        ])

        position = """&position
        %s
        /
        """ % pos

        s = tl.headtrim(general + interaction + cell + cutoff + position)
        tl.write(s, 'alm.in')
        write_vasp(
            'POSCAR-supercell', atoms, sort="True", direct=True, vasp5=True)

    def get_anphonin(self):
        m = self.m
        masses = m.getMassFromLabel(m.elements)

        general = """&general
        \tPREFIX = alm
        \tMODE = phonons
        \tFCSXML = alm.xml
        \tNKD = %s
        \tKD = %s
        \tMASS = %s
        /
        """ % (len(m.elements), m.toString(m.elements), m.toString(masses))

        cell = """&cell
        \t1.889726
        \t%s
        /
        """ % ('\n\t'.join([m.toString(m.atoms.cell[i]) for i in range(3)]))

        bp = m.bandpoints
        s = ""
        for i in range(len(m.bandpath) - 1):
            x1, x2 = m.bandpath[i], m.bandpath[i + 1]
            s += "  %s " % x1[0] + m.toString(bp[x1]) + \
                " %s " % x2[0] + m.toString(bp[x2]) + " 101\n"

        kpoint = """&kpoint
        \t1
        %s/
        """ % s
        tl.write(general + cell + kpoint, 'band.in')
        kpoint = """&kpoint
        \t2
        \t%s
        /
        """ % m.toString(m.kpoints)
        tl.write(general + cell + kpoint, 'dos.in')
        tl.write(general.replace('phonons', 'RTA') + cell + kpoint, 'tc.in')
