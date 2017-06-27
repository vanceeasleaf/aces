# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2017-06-18 22:24:41
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-26 22:28:03

from aces.tools import mkdir, cd,  mkcd, shell_exec,\
    write, pwd, ls, to_txt
import aces.config as config
from ase import io, Atoms
from aces.io.vasp import writevasp
from aces.runners import Runner
from aces.graph import plot
import numpy as np
from aces.runners.minimize import minimize as minimize_input
from aces.runners.phonopy import runner as PRunner
from importlib import import_module as im
import time
from numpy.linalg import norm
from aces.io.phonopy.fc import nomalizeFC, readfc2
from ase.transport import TransportCalculator

from aces.f import capacity


class runner(Runner):

    def creatmini(self, m):
        print('creatmini')
        m.home = pwd()
        assert m.home != ''
        mkdir('minimize')
        cd('minimize')
        minimize_input(m)
        write(time.strftime('%Y-%m-%d %H:%M:%S'), 'done')
        cd('..')
        if m.engine == "lammps":
            return m.dump2POSCAR(m.home + '/minimize/range')
        else:

            return io.read(m.home + '/minimize/CONTCAR')

    def test(self):
        dm = .1
        omega = np.arange(dm, 60, dm)  # THz
        factor = 1e12**2 * 1e-20 * 1e-3 / 1.6e-19 / 6.23e23
        energies = (omega * 2.0 * np.pi)**2 * factor
        # energies=np.arange(0,10,.01)
        h = -np.array((-2, 1, 0, 1, -2, 1, 0, 1, -2)).reshape((3, 3))
        h1 = -np.array((-2, 1, 1, -2)).reshape((2, 2))
        # x=1.0/np.sqrt(2)
        # h1=h=-np.array((-2,x,0,0,x,-1,x,0,0,x,-2,x,0,0,x,-1)).reshape((4,4))
        #energies = np.arange(-3, 3, 0.1)
        calc = TransportCalculator(h=h, h1=h1, energies=energies, dos=True)
        T = calc.get_transmission()
        # print T
        dos = calc.get_dos() * omega
        plot([omega, 'Frequency (THz)'], [T, 'Transmission'],
             'test_green_transmission.png')
        plot([omega, 'Frequency (THz)'], [
             dos, 'Phonon Density of State'], 'test_green_dos.png')
    """
	def collect(self):
		leadm=self.preLead()
		fclead=self.fc('lead')
		fccenter=self.fc('center')
		#write(np.around(fc[:,:,0,0],3),'orig_forces')
		n=leadm.hatom
		fccenter[:n,-n:]=0
		fccenter[-n:,:n]=0
		#write(np.around(fccenter[:,:,0,0],3),'fccenter')
		fclead=fclead[:2*n][:2*n]
		fccenter=self.reshape(fccenter)
		fclead=self.reshape(fclead)
		return fccenter,fclead
	"""

    def collect(self):

        fccenter, fcleft, fcright = self.fc('center')
        fccenter = self.reshape(fccenter)

        fcleft = self.reshape(fcleft)

        fcright = self.reshape(fcright)
        write(np.around(fcleft, 3), 'fcleft')
        write(np.around(fcright, 3), 'fcright')
        write(np.around(fccenter, 3), 'fccenter')
        # fcleft1=self.fc('leftlead')
        # fcleft1=self.reshape(fcleft1)
        # fcright1=self.fc('rightlead')
        # fcright1=self.reshape(fcright1)
        #
        # don't ingore any elements in the array to print
        # np.set_printoptions(threshold='nan')
        # write((fcleft1-fcleft)*1e6,'fcleft1')
        # write(fcright1-fcright,'fcright1')
        return fccenter, fcleft, fcright

    def gettrans(self):
        print("Reading in force constants...")
        # if not exists("fcbin.npz"):
        fccenter, fcleft, fcright = self.collect()
        np.savez(
            "fcbin.npz",
            fccenter=fccenter,
            fcleft=fcleft,
            fcright=fcright)
        print("Caching force constans")
        import os
        m = self.m
        os.system(config.mpirun + " " + str(m.nodes * \
                  m.procs) + " ae trans_cal >log.out")
        self.post()

    def reduce(self):
        files = ls("tmp/result.txt*")
        omega = []
        trans = []
        dos = []
        for file in files:
            print(file)
            # using this but loadtxt because there may be data that missing
            # some columns when there are nan
            result = self.read(file)
            omega = np.r_[omega, result[:, 0]]
            trans = np.r_[trans, result[:, 1]]
            dos = np.r_[dos, result[:, 2]]

        omega = np.array(omega).flatten().T
        f = omega.argsort()
        omega = omega[f]
        trans = np.array(trans).flatten().T[f]
        dos = np.array(dos).flatten().T[f]
        to_txt(['omega', 'trans', 'dos'], np.c_[
               omega, trans, dos], 'result.txt')

    def getcenter(self):
        centerm = self.preCenter()
        self.phonopy('center', centerm)

    def generate(self):
        self.m.xp = 1
        # leadm=self.preLead()
        # self.phonopy('lead',leadm)
        self.getcenter()
        self.run()

    def run(self):
        self.m.xp = 1
        self.getlead()
        self.gettrans()
        self.post()

    def read(self, file):
        return np.genfromtxt(file, delimiter='\t', skip_header=1)

    def post(self):
        # 1eV = 8049 cm^(-1) => 1000emV=8049 cm-1 => cm-1/meV=1000/8049
        # 1cm^(-1) = 3 * 10^(10) hz =>Hz*cm=1/3e10
        # a cm^-1=b THz =>a=b *1e12 Hz*cm
        # a meV = b cm^-1 => a = b cm-1/meV
        # omcm=omega*521.471？
        self.reduce()
        result = self.read('result.txt')
        omega = result[:, 0]
        trans = result[:, 1]
        dos = result[:, 2]
        omcm = omega * 1e12 * 1 / 3e10
        omme = omcm * 1e12 * 6.6260755e-34 / 1.6e-19 * 1000

        T = self.m.T
        centerm = self.preCenter()
        V = np.linalg.det(centerm.atoms.cell)
        c = capacity(omega, T, V)
        j = c * trans / 2.0 / np.pi
        dm = omega[1] - omega[0]
        kappa = j.cumsum() * dm
        to_txt(['Frequency (THz)',
                'Frequency (cm^-1)',
                'Frequency (meV)',
                'Phonon Transmission',
                'Phonon Density of State',
                'Mode Capacity (J/m^3/K)',
                'Mode Thermal Conductance (W/m^2/K)',
                'Accumulate Thermal Conductance (W/m^2/K)'],
               np.c_[omega,
                     omcm,
                     omme,
                     trans,
                     dos,
                     c,
                     j,
                     kappa],
               'transmission.txt')

        f = self.read('transmission.txt')
        #from aces.algorithm.smooth import savitzky_golay
        plot([f[:, 0], 'Frequency (THz)'], [
             f[:, 4], 'Phonon Density of State'], 'green_dos.png')
        plot([f[:, 0], 'Frequency (THz)'], [
             f[:, 3], 'Phonon Transmission'], 'green_transmission.png')
        #plot([f[:,0],'Frequency (THz)'],[savitzky_golay(f[:,3],11,3),'Phonon Transmission'],'smooth_transmission.png')
        plot([f[:, 0], 'Frequency (THz)'], [
             f[:, 6], 'Mode Thermal Conductance (W/m^2/K)'], 'green_mode_conductance.png')

    def reshape(self, fc):
        n, m = fc.shape[:2]
        fc = np.einsum('ikjl', fc).reshape([n * 3, m * 3])
        return fc

    def testfc(self):
        self.fc('center')

    def fc(self, dir):
        fc = readfc2(dir + '/FORCE_CONSTANTS')
        satoms = io.read(dir + '/SPOSCAR')

        # divide m_i*mj

        fc = nomalizeFC(fc, satoms)
        if not dir == 'center':
            atoms = io.read(dir + '/POSCAR')

            fc = self.rearangefc(fc, satoms, atoms)
        else:
            left = io.read(dir + '/POSCAR_left')
            right = io.read(dir + '/POSCAR_right')
            # why not record the order when generating POSCAR_left ? because
            # atoms order may be changed after write to POSCAR
            fc = self.rearangefc_center(fc, satoms, left, right)

        return fc

    def phonopy(self, dir, mm):
        # if exists(dir+'/FORCE_CONSTANTS'):
        #	return
        mkcd(dir)
        self.creatmini(mm)
        PRunner(mm).generate()
        cd('..')

    def preCenter(self):
        m = self.m
        import device.device as s
        leadm = self.preLead()
        u = s.Device(m, leadm, leadm)
        u.cores = m.cores
        u.__dict__ = dict(m.__dict__, **u.__dict__)
        return u

    def getEndFil(self, atoms, refpos, la, end="right", err=0.4):
        """get the index of the left or right end with depth=la

        [description]

        Arguments:
                atoms {[Atoms]} -- [description]
                refpos {[array[3]]} -- [the position of relative atom that is trait as 0 when counting for depth]
                la {[number]} -- [the depth]

        Keyword Arguments:
                end {str} -- [left or right] (default: {"right"})
                err {number} -- [the tolerance for atoms selectoin] (default: {0.4})

        Returns:
                [array[N]] -- [the filter of the ends]
        """
        ax = self.m.negfaxis
        errs = [0.0, 0, 0]
        errs[ax] = err
        offset = np.array([0.0, 0, 0])
        if end == "right":
            offset[ax] = la
            # to account for tilt cells
            uatoms = Atoms(
                'C',
                positions=[
                    np.array(
                        refpos -
                        offset) -
                    errs],
                cell=atoms.cell)
            fil = atoms.get_scaled_positions()[
                :, ax] > uatoms.get_scaled_positions()[
                0, ax]

        else:
            offset[ax] = -la
            uatoms = Atoms(
                'C',
                positions=[
                    np.array(
                        refpos -
                        offset) +
                    errs],
                cell=atoms.cell)
            fil = atoms.get_scaled_positions()[
                :, ax] < uatoms.get_scaled_positions()[
                0, ax]
        return fil, offset

    def findunit(self, atoms, refpos, la, end="right", err=0.4):
        """[find atoms within the period]

        [filter the atoms with x > rightx-la (if end=='right') and use la as lattice constant to build a new unitcell]

        Arguments:
                atoms {[Atoms]} -- [The center region atoms of NEGF, with scatter part and left and right lead part of 2 layers]
                refpos {[Array[3]]} -- [the position of right most atom]
                la {[Number]} -- [the lattice constant of result unitcell]

        Keyword Arguments:
                end {str} -- [which lead do you want ,can be 'left' and 'right'] (default: {"right"})
                err {number} -- [the tolerence to find an atom ] (default: {0.3})

        Returns:
                [Atoms] -- [The result unitcell]
        """

        fil, offset = self.getEndFil(atoms, refpos, la, end, err)
        patoms = atoms[fil]
        print(offset)
        # exclude the atoms that could be get though move one righer atom
        invalids = []
        for p in patoms:
            for i, q in enumerate(patoms):
                if norm(
                        q.position +
                        offset -
                        p.position) < err and q.symbol == p.symbol:
                    invalids.append(i)
        qatoms = Atoms()
        for i, q in enumerate(patoms):
            if i in invalids:
                continue
            qatoms.append(q)
        # check qatoms is a unitcell
        isunit = [False] * len(qatoms)
        for i, p in enumerate(qatoms):
            for q in atoms:
                if norm(
                        p.position -
                        offset -
                        q.position) < err and q.symbol == p.symbol:
                    isunit[i] = True
                    break
        print("is maching:", isunit)
        # the unitcell must exist
        #assert reduce(lambda a,b:a*b,isunit)
        return qatoms, np.array(isunit)

    def getlead(self):
        cd('center')
        left = self.findlead('left')
        right = self.findlead('right')
        cd('..')
        return
        from aces.io.vasp import writevasp
        mkcd('leftlead')
        writevasp(left, 'POSCAR')
        self.runlead()
        cd('..')
        mkcd('rightlead')
        writevasp(right, 'POSCAR')
        self.runlead()
        cd('..')

    def clear(self):
        shell_exec("rm fc*;rm tmp -r;rm *.png;rm transmission.txt result.txt")

    def findlead(self, end="left", err=0.3):
        """[from center POSCAR find periodic lead POSCAR]

        [the center POSCAR is expected to consist of scatter region and left lead region of 2 layers and
        right lead region of 2 layers. So we can find the lead POSCAR from the two ends with this algorithm.]

        Keyword Arguments:
                end {str} -- [which lead do you want ,can be 'left' and 'right'] (default: {"left"})
                err {number} -- [the tolerence to find an atom ] (default: {0.4})

        Raises:
                Exception -- [end is not set correctly]
        """

        atoms = io.read('SPOSCAR')
        # find the most right atom index
        if end == "right":
            last = -1
            next_last = -2
            next_last1 = -3
        elif end == "left":
            last = 0
            next_last = 1
            next_last1 = 2
        else:
            raise Exception("Unkown value of `end`")
        ax = self.m.negfaxis
        right_idx = atoms.positions[:, ax].argsort()[last]
        ratom = atoms[right_idx]
        print("%s end atom position:" % end, ratom.position)
        # the atoms that has similar y,z to ratom , in our situatoin there is
        # at least two such atoms.
        cross = ~(np.arange(3, dtype='int') == ax)
        fil = norm(
            atoms.positions[
                :,
                cross] - ratom.position[cross],
            axis=1) < err / 10
        # Don't forget the symbol must be the same.
        # must convert to array
        fil1 = np.array(atoms.get_chemical_symbols()) == ratom.symbol
        fil = fil * fil1
        atoms_inline = atoms[fil]
        assert len(atoms_inline) >= 2
        right_idx1 = atoms_inline.positions[:, ax].argsort()[next_last]
        # waring this index is of atoms_inline but not of atoms
        # catom is what we search
        catom = atoms_inline[right_idx1]
        print("next %s end similar atom position:" % end, catom.position)
        # the possible lattice constant of x direction
        la = np.abs(ratom.position[ax] - catom.position[ax])
        refpos = ratom.position
        unit, isunit = self.findunit(atoms, refpos, la, end, err)
        if (isunit == False).sum() > 1:
            next_last = next_last1
            right_idx1 = atoms_inline.positions[:, ax].argsort()[next_last]
            # waring this index is of atoms_inline but not of atoms
            # catom is what we search
            catom = atoms_inline[right_idx1]
            print(
                "next next %s end similar atom position:" %
                end, catom.position)
            la = np.abs(ratom.position[ax] - catom.position[ax])
            refpos = ratom.position
            unit, isunit = self.findunit(atoms, refpos, la, end, err)
            assert isunit.all()
        # the possible lattice constant of x direction
        rightx = refpos[ax]
        if end == "right":
            offset = [0, 0, 0]
            offset[ax] = -(rightx - la)
            unit.translate(offset)
        else:
            offset = [0, 0, 0]
            offset[ax] = -(rightx)
            unit.translate(offset)
        unit.cell = atoms.cell
        unit.cell[ax, ax] = la
        ys = unit.positions[:, (ax + 1) % 3]
        order = ys.argsort()
        unit = unit[order]
        from aces.io.vasp import writevasp
        writevasp(unit, 'POSCAR_%s' % end)
        return unit

    def runlead(self):
        """generate lead force constants

        TRICK!!
        > If hc1/hc2 are None, they are assumed to be identical to the coupling matrix elements between neareste neighbor principal layers in lead1/lead2.

        3 is important for the NEGF calculation ,if use 2, the regurlar fc and periodic fc is undistinguable,
        1 -1
        -1 1

        and the transimission is curve
        for example ,
        the fclead should be when lead layer=1,and two layer interaction is
        2 -1
        -1 2
        and we can get from 3 layer supercell by  fc[:2n,:2n] , this process will complete in rearangefc
        2 -1 -1
        -1 2 -1
        -1 -1 2

        """
        m = self.m
        s = im('aces.materials.graphene')
        mm = s.structure(dict(latx=1, laty=1, latz=1, xp=1, yp=0, zp=0))
        mm.atoms = io.read("POSCAR")

        mm.supercell = [3, 1, 1]
        mm.phofc = True
        mm.__dict__ = dict(m.__dict__, **mm.__dict__)
        PRunner(mm).run()

    def preLead(self):
        m = self.m
        s = im('aces.materials.%s' % m.leads)
        lat = m.leadlat
        mm = s.structure(
            dict(
                latx=lat[0],
                laty=lat[1],
                latz=lat[2],
                xp=1,
                yp=1,
                zp=1))
        mm.dimension = m.dimension
        import device.lead as s
        u = s.Lead(mm)
        # u.cores=m.cores
        u.__dict__ = dict(m.__dict__, **u.__dict__)
        return u

    def rearangefc_center(self, fc, satoms, left, right):
        ax = self.m.negfaxis
        right_idx = satoms.positions[:, ax].argsort()[-1]

        x_right = satoms[right_idx].position[ax]
        left_idx = satoms.positions[:, ax].argsort()[0]

        x_left = satoms[left_idx].position[ax]
        # don't put these after satoms=satoms[order]
        lpos = satoms[left_idx].position
        rpos = satoms[right_idx].position
        # find the index of left part in satoms
        atoms = left.copy()
        # remove the offset
        offset = [0, 0, 0]
        offset[ax] = x_left
        atoms.translate(offset)
        lidx = self.match_idx(atoms, satoms)
        lidx = list(lidx)
        atoms.translate(left.cell[ax])
        lidx1 = self.match_idx(atoms, satoms)
        lidx1 = list(lidx1)
        # find the index of right part in satoms
        atoms = right.copy()
        offset = [0, 0, 0]
        offset[ax] = x_right
        atoms.translate(np.array(offset) - right.cell[ax])

        ridx = self.match_idx(atoms, satoms)
        ridx = list(ridx)
        atoms.translate(-right.cell[ax])
        ridx1 = self.match_idx(atoms, satoms)
        ridx1 = list(ridx1)
        # find the center index
        cidx = []
        for i in range(len(satoms)):
            if i in lidx:
                continue
            if i in ridx:
                continue
            cidx.append(i)
        order = lidx + cidx + ridx
        print(order)
        fccenter = fc[order][:, order]
        oatoms = satoms
        satoms = satoms.copy()
        satoms = satoms[order]
        writevasp(satoms, "POSCAR_ordered")

        cut = self.m.negfcut
        if cut > 0:
            self.elimiPeriodic(satoms, fccenter, cut, left, right, lpos, rpos)
        else:

            # the default assumption is that cell length is larger than interaction length,
            # so for | x x | x x | x x | the interaction of cell 1 and cell 3 is only from periodic and the periodic interction only apparear in between 1 and 3
            # F_13 <==> periodic interaction so if we want to eliminate
            # periodic interaction we only have to set F_13=0

            n1 = len(lidx)
            n2 = len(ridx)
            fccenter[:n1, -n2:] = 0
            fccenter[-n2:, :n1] = 0

        # with this method there is no need to calculate fc for left or right
        # again, it's very important for abinitial calculation.
        order = lidx + lidx1
        fcleft = fc[order][:, order]
        # this is for the situation of 2x1x1 supercell when the cell size if
        # pretty large
        if len(cidx) == 0 and cut > 0:
            satoms = oatoms[order]
            self.elimiPeriodic(satoms, fcleft, cut, left, left, lpos, rpos)
        order = ridx1 + ridx
        fcright = fc[order][:, order]
        if len(cidx) == 0 and cut > 0:
            satoms = oatoms[order]
            self.elimiPeriodic(satoms, fcright, cut, right, right, lpos, rpos)
        return fccenter, fcleft, fcright

    def elimiPeriodic(self, satoms, fc, cut, left, right, lpos, rpos):
        """eliminate periodic effect ,which is very important

        # if the FORCE_CONSTANTS is obtain from 2x1x1 supercell it's also possible to get correct fc without periodic effects
        # for example if we have only nearest interaction in  | x x | x x | then the 1 and 4 atom interaction <==> periodic interaction
        # that is , periodic interaction only present between 1 and 4 ,and 1 and 4 atoms has not in-cell interaction.
        # negfcut is the max distance from ends to filter those only has periodic interaction but not in-cell insteraction
        # we have to garentee that there is not such atoms that interact both in-cell and periodic, that's physically wrong when getting fc.


        Arguments:
                satoms {[Atoms]} -- [super cell ]
                fc {[fc[N][N][3][3]]} -- [force constant ,N=len(satoms)]
                cut {[number]} -- [max interaction distance]
                left {[Atoms]} -- [left cell]
                right {[Atoms]} -- [right cell]
                la {[nubmer]} -- [description]
        """
        ax = self.m.negfaxis
        la = min(cut, norm(left.cell[ax]))
        lfil, offset = self.getEndFil(satoms, lpos, la, "left")
        la = min(cut, norm(right.cell[ax]))
        rfil, offset = self.getEndFil(satoms, rpos, la, "right")
        fc[lfil][:, rfil] = 0
        fc[rfil][:, lfil] = 0

    def match_idx(self, atoms, satoms, err=0.4):
        """find the indexes of atoms in satoms with the same position

        [description]

        Arguments:
                atoms {[type]} -- [description]
                satoms {[type]} -- [description]
                err {Number} -- tolerance for position compare
        """
        ax = self.m.negfaxis
        natom = len(atoms)
        idx = -np.ones(natom, dtype="int")
        for i, p in enumerate(atoms):
            for j, q in enumerate(satoms):
                if norm(
                        q.position -
                        p.position) < err and q.symbol == p.symbol:
                    idx[i] = j
                if norm(
                        q.position +
                        satoms.cell[ax] -
                        p.position) < err and q.symbol == p.symbol:
                    idx[i] = j
                    break
        # there must be a match atom
        assert (idx >= 0).all()

        return idx

    def rearangefc(self, fc, satoms, atoms):
        """[ reorder fc by satoms-> atoms +atoms +atoms]

        the order of atoms in SPOSCAR is not [atoms in left+ atoms in right] but may be random
        for example ,the atoms in POSCAR is C4N4 then the order in SPOSCAR is C8N8 however we need C4N4+C4N4

        if we get the order = [0,4,1,5,2,6,3,7] ,which means newfc[1]=fc[4]
        then we have newfc[i]=fc[order[i]]=fc[order][i] => newfc=fc[order]
        more acurately because fc.dimension=2 we must have newfc=fc[order][:,order]

        Arguments:
                fc {[type]} -- [description]
                stoms {[type]} -- [description]
                atoms {[type]} -- [description]

        Returns:
                [type] -- [description]
        """

        #from aces.f import mapatoms,writefc2
        # pos,order=mapatoms(atoms,old)
        # old.write('old.xyz')
        # atoms[order].write('new.xyz')

        # find the index of left part in satoms
        ax = self.m.negfaxis
        lidx = self.match_idx(atoms, satoms)
        lidx = list(lidx)

        # find the index of right part in satoms
        ratoms = atoms.copy()
        ratoms.translate(satoms.cell[ax] - atoms.cell[ax])
        ridx = self.match_idx(ratoms, satoms)
        ridx = list(ridx)
        # concat them
        cidx = []
        for i in range(len(satoms)):
            if i in lidx:
                continue
            if i in ridx:
                continue
            cidx.append(i)
        order = lidx + cidx + ridx
        print(order)
        print(satoms.positions[order])
        newfc = fc[order][:, order]
        # x[:2,:2]==x[:2][:,:2]
        # but x[[0,1],[0,1]]==[x[0,0],x[1,1]] !=x[[0,1]][:,[0,1]] that's why
        # the order and :2len is not the same
        return newfc[:2 * len(atoms), :2 * len(atoms)]
        # generated from aces.io.vasp.writevasp to record original order of atoms
        # order=np.loadtxt(dir+'/POSCARswap').astype(np.int)
        # writefc2(fc[order][:,order],'fc')
        # return fc[order][:,order]
