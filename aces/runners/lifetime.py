# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2017-06-18 22:47:37
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-18 22:46:30

from aces.runners import Runner
from aces.runners.correlation import runner as Crun
from aces.runners.phonopy import runner as Prun
from aces.graph import plot, imshow
import numpy as np
import aces.tools as tl
hbar = 6.6260755e-34 / 3.14159 / 2.0
kb = 1.3806488e-23

# Bose-Einstein Distribution


def BE(w, T):
    w = np.array(w)
    t = hbar * w / kb / T
    # return np.exp(-t)
    return np.nan_to_num(1.0 / (np.exp(t) - 1.0))


class runner(Runner):

    def postp(self):
        prun = Prun(self.m)
        prun.postp()

    def generate(self):
        prun = Prun(self.m)
        prun.run()
        self.fc()
        self.vfc()
        m = self.m
        m.usephana = False
        # mkdir('phi')
        # main=pwd()
        # for i in range(m.nseed):
        # mkcd('phi/%d'%i)
        # self.corr()
        # cd(main)
        self.corr()
        self.nma()
        # self.sed()
        # self.band()

    def corr(self):
        crun = Crun(self.m)
        crun.run()
        # crun.vd.life_yaml(correlation_supercell=self.m.correlation_supercell)
        # self.drawlifetime()

    def fc(self):
        c = self.m.correlation_supercell
        q = []
        u = [int(x / 2) * 2 + 1 for x in c]
        for i in range(u[0]):
            for j in range(u[1]):
                for k in range(u[2]):
                    b = np.array([
                        float(i - c[0] / 2) / c[0],
                        float(j - c[1] / 2) / c[1],
                        float(k - c[2] / 2) / c[2]
                    ])
                    q.append(b)
        m = self.m

        tl.mkcd('qpoints')
        tl.cp('../FORCE_CONSTANTS', '.')
        tl.cp('../disp.yaml', '.')
        tl.cp('../POSCAR', '.')
        Prun(m).getqpoints(q)
        tl.cd('..')

    def vfc(self):
        c = self.m.correlation_supercell
        q = []
        u = [int(x / 2) * 2 + 1 for x in c]
        for i in range(u[0]):
            for j in range(u[1]):
                for k in range(u[2]):
                    b = np.array([
                        float(i - c[0] / 2) / c[0],
                        float(j - c[1] / 2) / c[1],
                        float(k - c[2] / 2) / c[2]
                    ])
                    q.append(b)
        m = self.m

        tl.mkcd('vqpoints')
        tl.cp('../FORCE_CONSTANTS', '.')
        tl.cp('../disp.yaml', '.')
        tl.cp('../POSCAR', '.')
        Prun(m).getvqpoints(q)
        tl.cd('..')

    def dos(self):
        crun = Crun(self.m)
        crun.dos()

    def nmaq(self, k=[0, 0, 0], test=False):
        from aces.runners.vdos import vdos
        correlation_supercell = self.m.correlation_supercell
        vdos(self.m.timestep).lifenmaq(
            k=k, correlation_supercell=correlation_supercell, test=test)

    def nma(self):
        from aces.runners.vdos import vdos
        correlation_supercell = self.m.correlation_supercell
        vdos(self.m.timestep).lifenma(
            correlation_supercell=correlation_supercell)

    def allnma(self):
        from aces.runners.vdos import vdos
        vdos(self.m.timestep).allnma()

    def allnmaq(self, k=[0, 0, 0]):
        from aces.runners.vdos import vdos
        vdos(self.m.timestep).allnmaq(k=k)

    def sed(self):
        from aces.runners.vdos import vdos
        correlation_supercell = self.m.correlation_supercell
        vdos(self.m.timestep).lifesed(
            correlation_supercell=correlation_supercell)

    def band(self):
        from aces.runners.vdos import vdos
        correlation_supercell = self.m.correlation_supercell
        vdos(self.m.timestep).sed_band(
            correlation_supercell=correlation_supercell)
        self.drawband()

    def drawband(self):
        sed = np.load('sed.npy')
        # n = len(sed)
        # x = np.arange(n)
        w = np.load('wtick.npy')
        filter = w < 60
        sed = sed[:, filter]
        # X,W=np.meshgrid(x,w)
        # scatter(X,W,sed,'Wave Vector',
        # 'Frequency (THz)','sed.png',marker_size=2,marker='s')
        imshow(np.log(sed.T), 'sed.png', extent=[0, 1, 0, 1])

    def reducephi(self, nd, *pn):
        from aces.qpointsyaml import phononyaml
        pya = phononyaml("qpoints/qpoints.yaml")
        nq = pya.nqpoint
        nbr = pya.nbranch

        import h5py
        nma = h5py.File('lifenma.h5')
        q = np.array(nma['/q/0/0'])
        n = len(q)
        v = np.zeros([nq * nbr, n])

        def onedir(dir, num):
            print(dir)
            nma = h5py.File(dir)
            for i in range(nq):
                for j in range(nbr):
                    node = '/%d/%d' % (i, j)
                    v[i * nbr + j, :] += np.array(nma[node]) * num
            return num

        N = 0
        for i in range(nd):
            dir, num = pn[i * 2:i * 2 + 2]
            N += onedir(dir, num)
        v /= N
        phis = h5py.File('phis.h5')
        for i in range(nq):
            for j in range(nbr):
                node = '/%d/%d' % (i, j)
                phis[node] = v[i * nbr + j, :]

    def reduce(self, nd, *pn):
        from aces.qpointsyaml import phononyaml
        pya = phononyaml("qpoints/qpoints.yaml")
        nq = pya.nqpoint
        nbr = pya.nbranch
        import h5py
        nma = h5py.File('lifenma.h5')
        q = np.array(nma['/q/0/0'])
        n = len(q)
        v = np.zeros([nq * nbr, n])

        def onedir(dir, num):
            for idir in num:
                print(idir)
                nma = h5py.File(dir + '/%d/lifenma.h5' % idir)
                for i in range(nq):
                    for j in range(nbr):
                        node = '/q/%d/%d' % (i, j)
                        v[i * nbr + j, :] += np.array(nma[node])
            return len(num)

        N = 0
        for i in range(nd):
            dir, num = pn[i * 2:i * 2 + 2]
            N += onedir(dir, num)
        # N+=onedir('../q3.1',30)
        # N+=onedir('../q3.2',90)
        v /= N
        if tl.exists('phis.h5'):
            tl.passthru('rm phis.h5')
        phis = h5py.File('phis.h5')
        for i in range(nq):
            for j in range(nbr):
                node = '/%d/%d' % (i, j)
                phis[node] = v[i * nbr + j, :]

    def drawlifetime(self):
        a = np.loadtxt('allnma.txt', skiprows=1)
        om = a[:, 3]
        tao = a[:, 6]  # np.abs(1/a[:,5])
        tao[np.abs(om) < 1e-6] = 0.0
        n = len(tao)

        filter = tao < 1e5
        om = om[filter]
        tao = tao[filter] / 6.28
        plot(
            (om, 'Frequency (THz)'), (tao, 'Relaxation Time (ps)'),
            'tao_freq.png',
            grid=True,
            scatter=True,
            logy=True)
        tl.to_txt(['freq', 'tao'], np.c_[om, tao], 'tao_freq.txt')

        v = np.loadtxt('vqpoints/v.txt')[:n, 4]
        v = v[filter]
        l = v * tao
        tl.to_txt(['freq', 'lamda'], np.c_[om[om.argsort()], l[om.argsort()]],
                  'lamda_freq.txt')

        plot(
            (om, 'Frequency (THz)'), (l, 'Mean Free Path (Angstrom)'),
            'lamda_freq.png',
            grid=True,
            scatter=True,
            logy=True)
        T = self.m.T
        w = om * 1e12 * 2.0 * np.pi
        from ase import io
        atoms = io.read('POSCAR')
        cs = self.m.correlation_supercell

        V = np.linalg.det(atoms.cell * cs)
        m = self.m
        # print m.enforceThick,m.thick,atoms.cell[2,2]
        if m.enforceThick:
            V *= m.thick / atoms.cell[2, 2]
        c = hbar * w * (BE(w, T + 0.005) - BE(w, T - 0.005)) * 100.0 / V * 1e30
        tl.to_txt(['freq',
                   'capacity(J/K)'], np.c_[om[om.argsort()], c[om.argsort()]],
                  'capacity_freq.txt')
        plot(
            (om[om.argsort()], 'Frequency (THz)'),
            (c[om.argsort()], 'Mode Specific Heat (J/K)'),
            'capacity_freq.png',
            grid=True,
            linewidth=2)
        tl.to_txt(['freq', 'cumsumcapacity'],
                  np.c_[om[om.argsort()], c[om.argsort()].cumsum()],
                  'cumsumcapacity_freq.txt')
        plot(
            (om[om.argsort()], 'Frequency (THz)'),
            (c[om.argsort()].cumsum(), 'Acummulate Specific Heat (J/K)'),
            'cumsumcapacity_freq.png',
            grid=True,
            linewidth=2)

        k = l * 1e-10 * c * v * 1e-10 / 1e-12
        tl.to_txt(['freq', 'kappa'], np.c_[om[om.argsort()], k[om.argsort()]],
                  'kappa_freq.txt')
        plot(
            (om, 'Frequency (THz)'), (k, 'Mode Themal Conductivity (W/mK)'),
            'kappa_freq.png',
            grid=True,
            scatter=True,
            logy=True,
            logx=False)

        tl.to_txt(['freq', 'cumsumkappa'],
                  np.c_[om[om.argsort()], k[om.argsort()].cumsum()],
                  'cumsumkappa_freq.txt')
        plot(
            (om[om.argsort()], 'Frequency (THz)'),
            (k[om.argsort()].cumsum(),
             'Acummulate Themal Conductivity (W/mK)'),
            'cumsumkappa_freq.png',
            grid=True,
            linewidth=2)

        tl.to_txt(['lamda', 'cumsumkappa'],
                  np.c_[l[l.argsort()], k[l.argsort()].cumsum()],
                  'cumsumkappa_lamda.txt')
        plot(
            (l[l.argsort()], 'Mean Free Path (Angstrom)'),
            (k[l.argsort()].cumsum(), 'Acummulate Themal Conductivity (W/mK)'),
            'cumsumkappa_lamda.png',
            grid=True,
            linewidth=2)
