# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2017-06-16 20:09:09
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-27 16:02:34

from aces.tools import mkdir, mv, cd, cp, mkcd, shell_exec,\
    exists, write, passthru, toString, pwd, debug, ls, parseyaml
import aces.config as config
from aces.binary import pr
from aces.runners import Runner
from aces.graph import plot, series, pl, fig
from aces.script.vasprun import exe as lammpsvasprun
import aces.script.vasprun as vasprun
import time
import numpy as np
from aces.io.phonopy.bandplot import plotband, plotbanddos
from aces.io.phonopy.meshyaml import meshyaml
from aces.io.phonopy.fc import readfc2
from aces.pbs.jobManager import jobManager, th, pbs
from aces.io.vasp import writePOTCAR, writevasp, parseVasprun
from ase import io
from lxml import etree
from scanf import sscanf


class runner(Runner):

    def minimizePOSCAR(self):
        m = self.m
        if m.engine == "lammps":

            m.dump2POSCAR(m.home + '/minimize/range', rotate=True)

        elif m.engine == "vasp":
            cp(m.home + '/minimize/CONTCAR', 'POSCAR')

    def optimize(self):
        mkcd('optimize')
        cp('../minimize/POSCAR', '.')
        atoms = io.read('POSCAR')
        for i in range(100):
            dir = "%i" % i
            mkcd(dir)
            writevasp(atoms)
            forces, stress, energy = self.energyForce()
            pos = atoms.get_scaled_positions()
            pos += forces * 0.01

    def energyForce(self):
        self.getVaspRun_vasp()
        forces = parseVasprun('forces')
        stress = parseVasprun('stress')
        c = shell_exec("grep TOTEN OUTCAR|tail -1")
        energy = sscanf(c, "free  energy   TOTEN  =      %f eV")[0]
        return forces, stress, energy

    def cs(self):
        from aces.cs import runner
        runner(NAH=2).run()
        self.check('csfc2')

    def check1(self, filename='FORCE_CONSTANTS'):
        ref = io.read('SPOSCAR')
        fc2 = readfc2(filename)
        np.set_printoptions(precision=2, suppress=True)
        files = ['dir_POSCAR-001']
        vasprunxml = "dir_SPOSCAR/vasprun.xml"
        if exists(vasprunxml):
            vasprun = etree.iterparse(vasprunxml, tag='varray')
            forces0 = parseVasprun(vasprun, 'forces')
            print(forces0.max())
        else:
            forces0 = 0.0
        for file in files:
            print(file)
            POSCAR = 'dirs/%s/POSCAR' % file
            vasprunxml = "dirs/%s/vasprun.xml" % file
            atoms = io.read(POSCAR)
            u = atoms.positions - ref.positions
            f = -np.einsum('ijkl,jl', fc2, u)

            vasprun = etree.iterparse(vasprunxml, tag='varray')
            forces = parseVasprun(vasprun, 'forces') - forces0
            print(np.abs(f).max(), "\n")
            print(np.abs(forces - f).max())
            print(np.allclose(f, forces, atol=1e-2))

    def check(self, filename='FORCE_CONSTANTS'):
        ref = io.read('SPOSCAR')
        files = shell_exec("ls dirs").split('\n')
        fc2 = readfc2(filename)
        np.set_printoptions(precision=2, suppress=True)
        vasprunxml = "dir_SPOSCAR/vasprun.xml"
        if exists(vasprunxml):
            vasprun = etree.iterparse(vasprunxml, tag='varray')
            forces0 = parseVasprun(vasprun, 'forces')
            print(forces0.max())
        else:
            forces0 = 0.0
        for file in files:
            print(file)
            POSCAR = 'dirs/%s/POSCAR' % file
            vasprunxml = "dirs/%s/vasprun.xml" % file
            atoms = io.read(POSCAR)
            u = atoms.positions - ref.positions
            f = -np.einsum('ijkl,jl', fc2, u)

            vasprun = etree.iterparse(vasprunxml, tag='varray')
            forces = parseVasprun(vasprun, 'forces') - forces0
            print(np.abs(f).max(), "\n")
            print(np.abs(forces - f).max())
            print(np.allclose(f, forces, atol=1e-2))

    def stub(self):
        files = shell_exec("ls dirs").split('\n')
        files = map(lambda x: x.replace('dir_', ''), files)
        fc2 = readfc2('fc2')
        for file in files:
            ref = io.read('SPOSCAR')
            a = 'dirs/dir_' + str(file)
            atoms = io.read(a + "/POSCAR")
            u = atoms.positions - ref.positions
            f = -np.einsum('ijkl,jl', fc2, u)
            forces = ""
            for force in f:
                forces += "<v>  %f %f %f </v>\n" % tuple(force)
            vasprun = '<root><calculation><varray name="forces" >\n'
            vasprun += forces
            vasprun += '</varray></calculation></root>\n'
            write(vasprun, a + "/vasprun.xml")

    def force_constant(self, files):
        cmd = config.phonopy + "-f "
        if exists("dir_SPOSCAR/vasprun.xml"):
            cmd = config.phonopy + "--fz dir_SPOSCAR/vasprun.xml "
        for file in files:
            dir = "dirs/dir_" + file
            cmd += dir + '/vasprun.xml '
        # generate FORCE_SETS
        passthru(cmd)
        m = self.m
        # Create FORCE_CONSTANTS
        passthru(config.phonopy + "--tolerance=1e-4  --writefc --dim='%s'" %
                 (m.dim))

    def fc2(self):
        files = shell_exec("ls dirs").split('\n')
        files = map(lambda x: x.replace('dir_', ''), files)
        # when the number of files >1000, the order is wrong ,POSCAR-001,
        # POSCAR-1500 ,POSCAR-159
        files.sort(lambda x, y: int(x.split('-')[1]) - int(y.split('-')[1]))
        self.force_constant(files)

    def generate_meshconf(self):
        # generate mesh.conf
        m = self.m

        mesh = """DIM = %s
        ATOM_NAME = %s
        MP = %s
        EIGENVECTORS=.TRUE.
        FORCE_CONSTANTS = READ
        MESH_SYMMETRY = .FALSE.
        PRIMITIVE_AXIS = %s
        """ % (m.dim, ' '.join(m.elements), ' '.join(map(str, m.kpoints)),
               toString(m.premitive.flatten()))
        mesh = mesh.replace(r'^\s+', '')
        write(mesh, 'mesh.conf')

    def generate_vconf(self):
        # generate v.conf
        m = self.m

        mesh = """DIM = %s
        ATOM_NAME = %s
        MP = %s
        FORCE_CONSTANTS = READ
        MESH_SYMMETRY = .FALSE.
        GROUP_VELOCITY=.TRUE.
        PRIMITIVE_AXIS = %s
        """ % (m.dim, ' '.join(m.elements), ' '.join(map(str, m.kpoints)),
               toString(m.premitive.flatten()))
        mesh = mesh.replace(r'^\s+', '')
        write(mesh, 'v.conf')

    def generate_qconf(self, q):
        # generate q.conf
        m = self.m

        mesh = """DIM = %s
        ATOM_NAME = %s
        FORCE_CONSTANTS = READ
        EIGENVECTORS=.TRUE.
        QPOINTS=.TRUE.
        PRIMITIVE_AXIS = %s
        """ % (m.dim, ' '.join(m.elements), toString(m.premitive.flatten()))
        mesh = mesh.replace(r'^\s+', '')
        write(mesh, 'q.conf')
        s = "%s\n" % len(q)
        for qq in q:
            s += "%s\n" % toString(qq)
        write(s, 'QPOINTS')

    def generate_vqconf(self, q):
        # generate q.conf
        m = self.m

        mesh = """DIM = %s
        ATOM_NAME = %s
        FORCE_CONSTANTS = READ
        GROUP_VELOCITY=.TRUE.
        QPOINTS=.TRUE.
        PRIMITIVE_AXIS = %s
        """ % (m.dim, ' '.join(m.elements), toString(m.premitive.flatten()))
        mesh = mesh.replace(r'^\s+', '')
        write(mesh, 'q.conf')
        s = "%s\n" % len(q)
        for qq in q:
            s += "%s\n" % toString(qq)
        write(s, 'QPOINTS')

    def generate_supercells(self):
        m = self.m
        # generate supercells

        passthru(config.phonopy + "--tolerance=1e-4  -d --dim='%s'" % (m.dim))

    def writeINCAR(self):
        m = self.m
        npar = 1
        for i in range(1, int(np.sqrt(m.cores)) + 1):
            if m.cores % i == 0:
                npar = i
        if m.ispin:
            ispin = "ISPIN=2"
        else:
            ispin = ""
        if m.soc:
            soc = "LSORBIT=T"
        else:
            soc = ""
        if m.isym:
            sym = "ISYM = 1"
        else:
            sym = "ISYM = 0"
        s = """SYSTEM=calculate energy
        PREC = High
        IBRION = -1
        ENCUT = %f
        EDIFF = 1.0e-8
        ISMEAR = %d; SIGMA = 0.01
        IALGO = 38
        LREAL = .FALSE.
        ADDGRID = .TRUE.
        LWAVE = .FALSE.
        LCHARG = .FALSE.
        NPAR = %d
        %s
        %s
        %s
        """ % (self.m.ecut, m.ismear, npar, sym, ispin, soc)
        if m.vdw:
            s += """\nIVDW = 1
            VDW_RADIUS = 50
            VDW_S6 = 0.75
            VDW_SR = 1.00
            VDW_SCALING = 0.75
            VDW_D = 20.0
            VDW_C6 = 63.540 31.50
            VDW_R0 = 1.898 1.892
            """
        s = s.replace(r'^\s+', '')
        write(s, 'INCAR')

    def getVaspRun_vasp(self):

        self.writeINCAR()
        m = self.m
        writePOTCAR(m, m.elements)

        if (m.kpointspath):
            cp(m.kpointspath, "KPOINTS")
        else:
            from aces.io.vasp import writeKPOINTS
            writeKPOINTS(m.ekpoints)
        if 'jm' in self.__dict__:
            if not m.th:
                path = pwd()
                if m.queue == "q3.4":
                    pb = pbs(
                        queue=m.queue,
                        nodes=12,
                        procs=1,
                        disp=m.pbsname,
                        path=path,
                        content=config.mpirun + " 12 " + config.vasp +
                        ' >log.out')
                else:
                    pb = pbs(
                        queue=m.queue,
                        nodes=1,
                        procs=12,
                        disp=m.pbsname,
                        path=path,
                        content=config.mpirun + " 12 " + config.vasp +
                        ' >log.out')
            else:
                path = pwd()
                pb = th(disp=m.pbsname, path=path)
            self.jm.reg(pb)

        else:
            shell_exec(config.mpirun + " %s " % m.cores + config.vasp +
                       ' >log.out')

    def getVaspRun_lammps(self):
        m = self.m
        if 'jm' in self.__dict__:
            path = pwd()
            pb = pbs(
                queue=m.queue,
                nodes=1,
                procs=4,
                disp=m.pbsname,
                path=path,
                content=config.python + vasprun.__file__ + ' >log.out')
            self.jm.reg(pb)
        else:
            shell_exec(config.python + vasprun.__file__ + ' >log.out')

    def thcode(self, files, put):
        s = ""
        for file in files:
            dir = "dirs/dir_" + file
            s += "cd %s\n" % (dir)
            s += "yhbatch -N 1 aces.pbs\n"
            s += "cd ../../\n"
        write(s, put + "/runall.sh")

    def getvasprun(self, files):
        m = self.m
        maindir = pwd()
        if m.engine == "vasp":
            calculator = self.getVaspRun_vasp
        elif m.engine == "lammps":
            calculator = self.getVaspRun_lammps
        self.jm = jobManager()
        for file in files:
            print(file)
            dir = "dirs/dir_" + file
            mkdir(dir)
            mv(file, dir + '/POSCAR')
            cd(dir)
            calculator()
            cd(maindir)
        self.jm.run()
        if m.th:
            mkdir(m.pbsname)
            self.thcode(files, m.pbsname)
            cp("dirs", m.pbsname)
            passthru("tar zcf %s.tar.gz %s" % (m.pbsname, m.pbsname))
        print('start check')
        self.jm.check()
        if m.engine == "lammps1":
            from multiprocessing.dummy import Pool
            pool = Pool()
            pool.map_async(lammpsvasprun, files)
            pool.close()
            pool.join()

    def runSPOSCAR(self):
        m = self.m
        maindir = pwd()
        file = "SPOSCAR"
        dir = "dir_" + file
        mkdir(dir)
        cp(file, dir + '/POSCAR')
        cd(dir)
        if m.engine == "vasp":
            self.getVaspRun_vasp()
        if m.engine == "lammps":
            self.getVaspRun_lammps()
        cd(maindir)

    def checkMinimize(self):
        import yaml
        data = yaml.load(open('disp.yaml').read())
        disps = [map(float, a['direction']) for a in data['displacements']]
        maindir = pwd()
        dirs = ls('dirs/dir_*')
        ii = 0
        L = np.linalg.norm
        # d,p,d1,p1=self.m.rot
        out = open('ccos.txt', 'w')
        for dir in dirs:
            cd(dir)
            f = open('dump.force')
            for i in range(9):
                f.next()
            for b in range(ii):
                f.next()
            line = f.next()
            line = line.split()
            force = np.array(map(float, line[1:4]))
            # force=RotateVector(force,d1,-p1)
            # force=RotateVector(force,d,-p)
            d = disps[i]
            ccos = force.dot(d) / L(force) / L(d)
            ii += 1
            print >> out, "%d\t%f" % (ii, ccos)
            cd(maindir)

    def run(self):
        m = self.m
        a = time.time()
        self.generate_supercells()
        debug('generate_supercells:%f s' % (time.time() - a))
        files = shell_exec("ls *-*").split('\n')
        assert len(files) > 0 and not files[0] == ""
        # self.runSPOSCAR()
        a = time.time()
        self.getvasprun(files)
        debug('getvasprun:%f s' % (time.time() - a))
        a = time.time()
        self.fc2()
        debug('force_constant:%f s' % (time.time() - a))

        if m.phofc:
            return self
        self.postp()

    def generate(self):

        self.minimizePOSCAR()
        self.run()

    def get_force_sets(self):
        files = shell_exec("ls dirs").split('\n')
        files = map(lambda x: x.replace('dir_', ''), files)
        self.force_constant(files)

    def postp(self):
        m = self.m
        if m.gamma_only:
            self.getDos()
            return
        self.getband()
        self.getDos()

        self.getbanddos()
        self.drawpr()
        self.getV()

    def getqpoints(self, q):
        self.generate_qconf(q)
        passthru(config.phonopy + "--tolerance=1e-4  q.conf")

    def getvqpoints(self, q):
        self.generate_vqconf(q)
        passthru(config.phonopy + "--tolerance=1e-4  q.conf")
        data = parseyaml('qpoints.yaml')
        file = open("v.txt", 'w')
        for phonon in data['phonon']:
            qp = phonon['q-position']
            for band in phonon['band']:
                frequency = band['frequency']
                v = np.array(band['group_velocity'])
                v = np.linalg.norm(v)
                print >> file, "%s\t%f\t%f" % ('\t'.join(map(str, qp)),
                                               frequency, v)
        file.close()
        v = np.loadtxt('v.txt')
        plot(
            (v[:, 3], 'Frequency (THz)'), (v[:, 4],
                                           'Group Velocity (Angstrom/ps)'),
            'v_freq.png',
            grid=True,
            scatter=True)

    def getDos(self):
        self.generate_meshconf()
        passthru(config.phonopy + "--tolerance=1e-4  --dos  mesh.conf")
        self.drawDos()

    def getV(self):
        if not exists('groupv'):
            mkdir('groupv')
        cd('groupv')
        cp('../FORCE_CONSTANTS', '.')
        cp('../POSCAR', '.')
        cp('../disp.yaml', '.')
        self.generate_vconf()
        passthru(config.phonopy + "--tolerance=1e-4    v.conf")
        self.drawV()
        cd('..')

    def drawV(self):
        data = parseyaml('mesh.yaml')
        file = open("v.txt", 'w')
        for phonon in data['phonon']:
            qp = phonon['q-position']
            for band in phonon['band']:
                frequency = band['frequency']
                v = np.array(band['group_velocity'])
                v = np.linalg.norm(v)
                print >> file, "%s\t%f\t%f" % ('\t'.join(map(str, qp)),
                                               frequency, v)
        file.close()
        v = np.loadtxt('v.txt')
        plot(
            (v[:, 3], 'Frequency (THz)'), (v[:, 4],
                                           'Group Velocity (Angstrom/ps)'),
            'v_freq.png',
            grid=True,
            scatter=True)

    def getband(self):
        self.generate_bandconf()
        passthru(config.phonopy + "--tolerance=1e-4  -s  band.conf")
        plotband(labels=' '.join(self.m.bandpath))

    def getbanddos(self):
        freq, pdos = self.getpdos()
        plotbanddos(
            freq=freq,
            dos=np.sum(pdos, axis=1),
            labels=' '.join(self.m.bandpath))

    def modulation(self):
        m = self.m
        conf = """
        DIM = %s
        MODULATION = 1 1 1, 0 0 0 0 1 0
        ATOM_NAME = %s
        FORCE_CONSTANTS = READ
        """ % (m.dim, ' '.join(m.elements))
        write(conf, 'modulation.conf')

        passthru(config.phonopy + "--tolerance=1e-4    modulation.conf")

    def animate(self):
        m = self.m
        conf = """
        DIM = %s
        ANIME = 0 5 20
        ANIME_TYPE = xyz
        ATOM_NAME = %s
        FORCE_CONSTANTS = READ
        """ % (m.dim, ' '.join(m.elements))
        write(conf, 'animate.conf')

        passthru(config.phonopy + "--tolerance=1e-4    animate.conf")

    def generate_bandconf(self):
        # generate mesh.conf
        m = self.m

        bp = m.bandpoints
        bpath = ' '.join([toString(bp[x]) for x in m.bandpath])

        band = """DIM = %s
        ATOM_NAME = %s
        BAND = %s
        BAND_POINTS = 101
        FORCE_CONSTANTS = READ
        PRIMITIVE_AXIS = %s
        """ % (m.dim, ' '.join(m.elements),
               bpath, toString(m.premitive.flatten()))
        band = band.replace(r'^\s+', '')
        write(band, 'band.conf')

    def getpdos(self):
        xx = np.loadtxt('partial_dos.dat', skiprows=1)
        freq = xx[:, 0]
        pdos = xx[:, 1:]
        return freq, pdos

    def drawDos(self):
        freq, pdos = self.getpdos()
        datas = [(freq, p, '') for p in pdos.T]
        series(
            'Frequency (THz)',
            'Partial Density of States',
            datas=datas,
            filename='partial_dos.png',
            legend=False,
            grid=True)

        plot(
            (freq, 'Frequency (THz)'), (np.sum(pdos, axis=1),
                                        'Density of States'),
            filename='total_dos.png')
        # calculate paticipation ratio

    def mesh(self):
        """ save mesh.yaml to mesh.npz

        [description]
        """
        data = meshyaml('mesh.yaml')
        np.savez('mesh', **data)

    def drawpr(self):
        pr()
        # plot
        xs = []
        ys = []
        for line in open('pr.txt'):
            x, y = map(float, line.split())
            xs.append(x)
            ys.append(y)
        write("%s" % (sum(ys) / len(ys)), "ave_pr.txt")
        with fig('Paticipation_ratio.png'):
            pl.plot(xs, ys, '.', color='red')
            pl.ylim([0.0, 1.0])
            pl.xlabel('Frequency (THz)')
            pl.ylabel('Paticipation Ratio')
