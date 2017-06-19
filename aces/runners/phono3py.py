# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2015-10-19 14:58:54
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-19 12:34:14

import aces.tools as tl
import aces.config as config
from ase.io import read
from aces.graph import plot
import numpy as np
from aces.runners.phonopy import runner as Runner


class runner(Runner):

    def force_constant(self, files):
        cmd = config.phono3py + "--cf3 "
        for file in files:
            dir = "dirs/dir_" + file
            cmd += dir + '/vasprun.xml '
        tl.write(cmd, 'post.sh')
        # generate FORCE_SETS
        tl.passthru("sh post.sh")
        m = self.m
        print("Create fc2.hdf and fc3.hdf")
        tl.passthru(config.phono3py + "-c POSCAR --dim='%s'" % (m.dim))

    def generate_supercells(self):
        m = self.m
        # generate supercells
        tl.passthru(config.phono3py + "-d --dim='%s' -c POSCAR" % (m.dim))

    def generate(self):
        m = self.m
        self.minimizePOSCAR()
        self.generate_supercells()
        files = tl.shell_exec("ls *-*").split('\n')
        self.getvasprun(files)
        self.force_constant(files)
        if m.pho3bte:
            self.sheng()
        else:
            self.caltc()
            self.dumpkappa()
        # Thermal conductivity calculation

    def caltc(self):
        m = self.m
        tl.passthru(config.phono3py + '  --dim="' + m.dim +
                    '" -c POSCAR --mesh="' + ' '.join(map(str, m.kpoints)) +
                    '" --fc3 --fc2  --thm --br')
        filename = "kappa-m%s.hdf5" % ''.join(map(str, m.kpoints))
        tl.passthru(config.pypath + 'kaccum   --mesh="' +
                    ' '.join(map(str, m.kpoints)) + '"  POSCAR ' + filename +
                    ' |tee kaccum.dat')

    def sheng(self):

        self.writeFC()
        m = self.m
        from aces.runners.shengbte import runner as shengbte
        a = shengbte(m)
        tl.mkcd('SHENG')
        tl.cp('../FORCE_CONSTANTS_2ND', '.')
        tl.cp('../FORCE_CONSTANTS_3RD', '.')
        a.getControl()
        # Thermal conductivity calculation
        print("START SHENGBTE...")
        # passthru(config.mpirun+" %s "%m.cores+config.shengbte)
        tl.passthru(config.shengbte)

    def writeFC(self):
        print("writing text FORCE CONSTANTS 2 from hdf5")
        import h5py
        f = h5py.File('fc2.hdf5')
        fc2 = f['fc2']
        from aces.f import disp2atoms, writefc2, writefc3
        writefc2(fc2, 'FORCE_CONSTANTS_2ND')
        f.close()
        print("writing text FORCE CONSTANTS 3 from hdf5")
        f = h5py.File('fc3.hdf5')
        # phono3py phi is bac , sheng bte phi is bca and shengbte thirdorder
        # need 100.0 scale
        fc3 = np.einsum(f['fc3'], [0, 2, 1, 3, 5, 4])
        atoms = read('POSCAR')
        satoms = disp2atoms('disp_fc3.yaml')
        writefc3(fc3, atoms, satoms, 'FORCE_CONSTANTS_3RD')

    def dumpkappa(self):
        import h5py
        m = self.m
        filename = "kappa-m%s.hdf5" % ''.join(map(str, m.kpoints))
        f = h5py.File(filename)
        kappa = f['kappa']
        T = f['temperature']
        tl.to_txt(['temperature', 'kappa'], np.c_[T, kappa[:, 0]], 'kappa.txt')
        plot(
            (np.array(T), 'Temperature (K)'),
            (np.array(kappa),
             'Thermal Conductivity (W/mK)'),
            filename='kT.png')
