# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2017-06-18 23:16:05
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-18 23:16:14

from aces.runners import Runner
from aces.runners.correlation import runner as Crun
from aces.runners.phonopy import runner as Prun
import aces.tools as tl


class runner(Runner):

    def corr(self):
        crun = Crun(self.m)
        crun.run()

    def generate(self):

        print("running phonopy to generate FORCE_SETS")
        prun = Prun(self.m)
        prun.run()
        print("running lammps to get dynaphopy.trajectory")
        self.corr()
        self.post()

    def post(self):
        self.dyinput()
        self.v2h5()
        self.vdos()

    def v2h5(self):
        if not tl.exists('dy.h5'):
            self.pix('-sv dy.h5')

    def dyinput(self):
        m = self.m
        f = open('in.dy', 'w')
        print >> f, "STRUCTURE FILE POSCAR"
        print >> f, tl.pwd() + "/POSCAR"
        print >> f, "FORCE SETS"
        print >> f, tl.pwd() + "/FORCE_SETS"
        # print >>f,"STRUCTURE FILE OUTCAR"
        # print >>f,pwd()+"/dynaphopy.lammpstrj"
        print >> f, "SUPERCELL MATRIX PHONOPY"
        print >> f, " %d 0 0 \n 0 %d 0\n 0 0 %d" % tuple(m.supercell)
        print >> f, "BANDS"
        bp = m.bandpath
        bpp = m.bandpoints
        for i in range(len(bp) - 1):
            u = "%s %s" % (tl.toString(bpp[bp[i]], ','),
                           tl.toString(bpp[bp[i + 1]], ','))
            print >> f, u
        f.close()

    def sed(self, k=[0, 0, 0]):
        self.pix("-q %s  -sw 'dvsed%s.txt' -r 0 50 2000" % (tl.toString(k),
                                                            str(k)))

    def vdos(self):
        self.pix("-sd dyvdos.txt")

    def nma(self, k=[0, 0, 0]):
        self.pix("-q %s  -sp 'dvnma%s.txt' -r 0 16 2000" % (tl.toString(k),
                                                            str(k)))

    def pix(self, cmd):
        if tl.exists('dy.h5'):
            v = "-lv dy.h5"
        else:
            v = "dynaphopy.lammpstrj"
        tl.passthru("dynaphopy in.dy %s -ts %f %s" % (v, self.m.timestep, cmd))

    def atomdisp(self, direction=[1, 0, 0]):
        self.pix(" -sad %s dydisplacement.txt" % (tl.toString(direction)))

    def life(self, k=[0, 0, 0]):
        self.pix("-pa -q %s --silent > dylife%s.txt" % tl.toString(k), str(k))

    def normalforce(self):
        self.pix("-sfc FORCE_CONSTANTS_DY")
