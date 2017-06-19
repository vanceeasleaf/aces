# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2016-12-24 18:47:10
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-19 12:36:00

from aces import config
import aces.tools as tl
from aces.runners.phonopy import runner as Runner
import time


class runner(Runner):

    def force_constant(self):
        cmd = config.phonopy + "--fc dirs/dir_POSCAR001/vasprun.xml"
        tl.passthru(cmd)

    def getVaspRun_vasp(self):
        m = self.m

        if m.isym:
            sym = "ISYM = 1"
        else:
            sym = "ISYM = 0"
        if m.ispin:
            spin = "ISPIN =2"
        else:
            spin = ""
        s = """SYSTEM=calculate energy
        PREC = High
        IBRION = 8
        ENCUT = %f
        EDIFF = 1.0e-8
        ISMEAR = %d; SIGMA = 0.01
        IALGO = 38
        LREAL = .FALSE.
        ADDGRID = .TRUE.
        LWAVE = .FALSE.
        LCHARG = .FALSE.
        %s
        %s
        """ % (self.m.ecut, m.ismear, sym, spin)
        s = s.replace(r'^\s+', '')
        tl.write(s, 'INCAR')
        m = self.m
        m.writePOTCAR()
        s = """A
        0
        Monkhorst-Pack
        %s
        0  0  0
        """ % ' '.join(map(str, m.ekpoints))
        s = s.replace(r'^\s+', '')
        tl.write(s, 'KPOINTS')
        if 'jm' in self.__dict__:
            from aces.jobManager import pbs
            path = tl.pwd()
            if m.queue == "q3.4":
                pb = pbs(
                    queue=m.queue,
                    nodes=12,
                    procs=4,
                    disp=m.pbsname,
                    path=path,
                    content=config.mpirun +
                    " 48 " +
                    config.vasp +
                    ' >log.out')
            else:
                pb = pbs(
                    queue=m.queue,
                    nodes=4,
                    procs=12,
                    disp=m.pbsname,
                    path=path,
                    content=config.mpirun +
                    " 48 " +
                    config.vasp +
                    ' >log.out')
            self.jm.reg(pb)
        else:
            tl.shell_exec(
                config.mpirun +
                " %s " %
                m.cores +
                config.vasp +
                ' >log.out')

    def generate(self):
        m = self.m
        self.minimizePOSCAR()
        a = time.time()
        self.generate_supercells()
        tl.debug('generate_supercells:%f s' % (time.time() - a))
        tl.shell_exec("rm *-*")
        tl.cp("SPOSCAR", "POSCAR001")
        a = time.time()
        files = ['POSCAR001']
        self.getvasprun(files)
        tl.debug('getvasprun:%f s' % (time.time() - a))
        a = time.time()
        self.force_constant()
        tl.debug('force_constant:%f s' % (time.time() - a))

        if m.phofc:
            return self
        self.postp()
