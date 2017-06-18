# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2015-10-20 22:19:29
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-18 22:40:43

import aces.tools as tl
import aces.config as config
from ase.io import read
from aces.runners import Runner
from aces.graph import plot
import numpy as np


class runner(Runner):

    def runcmd(self):
        exe = config.phonts + "  >log.out"
        return config.mpirun + " %s " % self.m.cores + exe

    def generate(self):
        m = self.m
        coordination = self.phontsAtoms()
        content0 = "species %d\n" % (len(m.elements)) + m.phontsmasses + """
D3_cutoff %f
kpoints %s 1
delta 0.005
numerical_2der T
numerical_3der T
iter_steps 3
pdos 0. 70. 200 10.
temp 60. 400. 10
AbInitio  T F
FP_interface LAMMPS
#phonons_only T
Lattice  1.0
%s
end
""" % (m.shengcut, m.toString(m.kpoints), coordination)
        tl.write(content0, 'phonons_input.dat')
        tl.passthru(config.phonts)  # generate many displacement files
        tl.mkdir('lammps')
        tl.cd('lammps')
        content = "units %s\n" % m.units
        content += """atom_style      charge
dimension       3
boundary        p p p
read_data       GENERIC
%s
%s
neighbor        1.1 bin
neigh_modify    every 1 delay 1 check yes
dump 1 all custom 1 *.dump id  fx fy fz
dump_modify 1 format "%%d %%30.20f %%30.20f %%30.20f"
dump_modify  1 sort id
run 0
""" % (m.masses, m.potential)
        tl.shell_exec("mv ../*.str .")
        strs = tl.shell_exec("ls *.str").split("\n")
        for str in strs:
            dir = str.replace("str", "dir")
            tl.mkdir(dir)
            tl.write(content.replace("GENERIC", str), dir + "/in")
            tl.mv(str, "%s/%s" % (dir, str))
            tl.cd(dir)
            tl.passthru(config.lammps + " <in >out.dat")
            tl.cd('..')
        dirs = tl.shell_exec("ls |grep dir").split("\n")
        for dir in dirs:
            print(dir)
            tl.cp(dir + "/0.dump", "../" + dir.replace("dir", "out"))
        tl.cp("1.0000.dir/out.dat", "../1.0000.out")
        tl.cd('..')
        content0 = content0.replace('AbInitio  T F', 'AbInitio  F T')
        tl.write(content0, 'phonons_input.dat')
        tl.passthru(config.phonts)

    def phontsAtoms(self):
        m = self.m
        m.dump2POSCAR('minimize/range')
        atoms = read('POSCAR')
        cell = atoms.get_cell()
        if not np.allclose(np.diag(np.diag(cell)), cell):
            raise Exception('phonts needs cell to be orthorgnal')
        content = "cell %f %f %f\n" % tuple(np.diag(cell))
        content += "natoms %d\n" % (len(atoms))
        content += "fractional\n"
        pos = atoms.get_scaled_positions()
        for i, atom in enumerate(atoms):
            content += "%s %s\n" % (atom.symbol, m.toString(pos[i]))
        return content

    def post(self):
        a = np.loadtxt('freq.dat')
        ks = a[:, 1:4]
        omega = a[:, 4:]
        b = np.loadtxt('phon_lifetime.dat')
        tao = b[:len(ks), 4:]

        v = np.loadtxt(open('group_vel.dat'))[:, 4:]
        n, m = v.shape
        v = v.reshape([n, m / 3, 3])
        v = np.linalg.norm(v, axis=2)
        plot(
            (omega.flatten(), 'Frequency (THz)'), (v.flatten(),
                                                   'Group Velocity (nm/ps)'),
            'v_freq.png',
            grid=True,
            scatter=True)
        tl.to_txt(['freq', 'vg'],
                  np.c_[omega.flatten(), v.flatten()], 'v_freq.txt')
        plot(
            (omega.flatten(), 'Frequency (THz)'), (tao.flatten(),
                                                   'Relaxation Time (ps)'),
            'tao_freq.png',
            grid=True,
            scatter=True,
            logy=True)
