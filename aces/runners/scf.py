# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2016-09-05 19:22:06
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-18 22:36:31

import aces.tools as tl
from aces.runners.phonopy import runner as Runner


class runner(Runner):

    def generate(self):
        tl.cp('minimize/POSCAR', '.')
        self.getVaspRun_vasp()

    def q(self):
        a = tl.shell_exec(
            "grep TOTEN OUTCAR |tail -1").split("=")[1].strip().replace("eV", "")
        print(self.m.ecut, a)
