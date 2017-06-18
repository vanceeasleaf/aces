# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2015-11-17 14:24:57
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-18 22:34:28


import aces.config as config
import aces.tools as tl
from aces.runners import Runner
from aces.lammpsdata import lammpsdata


class runner(Runner):

    def get_structure(self):
        atoms = self.m.atoms_from_dump('minimize/range')
        # atoms=atoms.repeat(self.m.correlation_supercell)
        a = lammpsdata(atoms, self.m.elements)
        a.writedata('strain_structure')

    def generate(self):
        m = self.m
        self.get_structure()
        f = open("strain.lmp", "w")
        print >>f, "units %s" % m.units
        print >>f, "dimension 3"
        pbcx = "f"
        if m.vStrain:
            pbcx = "s"
        pbcy = pbcz = 's'
        if m.xp == 1:
            pbcx = 'p'
        if m.yp == 1:
            pbcy = 'p'
        if m.zp == 1:
            pbcz = 'p'
        print >>f, "boundary %s %s %s" % (pbcx, pbcy, pbcz)
        print >>f, "atom_style atomic"
        print >>f, "read_data   strain_structure"
        print >>f, "change_box	all	boundary %s %s %s" % (pbcx, pbcy, pbcz)
        print >>f, "lattice fcc 5"  # needed to define the regions
        print >>f, "thermo %d" % m.dumpRate
        print >>f, "thermo_modify     lost warn"
        print >>f, m.masses
        print >>f, m.potential
        print >>f, "timestep %f" % m.timestep
        print >>f, "reset_timestep 0"

        box = m.box
        deta = m.deta
        wfix = m.wfix
        xlo, xhi, ylo, yhi, zlo, zhi, lx, ly, lz = box
        fixl1 = xlo - deta
        fixl2 = fixl1 + deta * wfix
        fixr2 = xhi + deta
        fixr1 = fixr2 - deta * wfix
        if m.vStrain:
            runTime, content = self.vDeform()
        else:
            runTime, content = self.deform()
        print >>f, "region	stayl	block   %s  %s INF  INF INF  INF units box" % (
            fixl1, fixl2)
        print >>f, "region	stayr	block   %s  %s INF INF   INF  INF units box" % (
            fixr1, fixr2)
        print >>f, "region   stay    union  2 stayl stayr"
        print >>f, "region	main	block   %s  %s INF INF   INF  INF units box" % (
            fixl2, fixr1)
        print >>f, "group   stayl    region  stayl"
        print >>f, "group   stayr    region  stayr"
        print >>f, "group   stay    region  stay"
        print >>f, "group   main    region  main"
        print >>f, "velocity stay set 0 0 0"
        print >>f, "velocity main create %f %d mom yes rot yes dist gaussian" % (
            m.T, m.seed)
        #print >>f,"velocity stay set NULL 0 0"
        #print >>f,"fix force stay setforce  NULL 0 0"
        print >>f, "fix getEqu  main  nvt temp %f %f %f" % (m.T, m.T, m.dtime)
        print >>f, "dump dump1 all atom %d dump.lammpstrj" % (
            max(runTime / 1000, 1))
        print >>f, "dump_modify  dump1 sort id"
        print >>f, "run %d" % m.equTime
        print >>f, "unfix getEqu"

        print >>f, "reset_timestep 0"
        print >>f, "fix nve main nve"
        print >>f, "compute    disp all displace/atom"
        print >>f, "compute    s1 all stress/atom NULL"
        print >>f, "compute    rr stayr com"
        print >>f, "compute    rl stayl com"
        print >>f, "fix  2  main temp/berendsen %f %f 1" % (m.T, m.T)
        print >>f, "fix s all ave/atom 1 %d %d c_disp[1] c_s1[1] c_s1[2] c_s1[3] c_s1[4] c_s1[5] c_s1[6]" % (
            m.strainStep, m.strainStep)
        print >>f, "dump 1 all custom %d dump.tensile id  type xs ys zs f_s[1] f_s[1] f_s[2] f_s[3] f_s[4] f_s[5] f_s[6] f_s[7]" % m.strainStep
        print >>f, "dump_modify  1 sort id"
        print >>f, "variable lx equal c_rr[1]-c_rl[1]"
        print >>f, "variable pxx equal pxx"
        print >>f, "fix p all ave/time 1 %d %d v_lx v_pxx file  strain_stress.txt" % (
            m.strainStep, m.strainStep)
        #print >>f,"dump lala main custom %s velocity.txt id type vx vy vz"%m.Cinterval
        print >>f, content
        f.close()
        tl.passthru(
            config.mpirun +
            "  %s " %
            self.m.cores +
            config.lammps +
            " <strain.lmp  >out.dat")
        self.post()

    def deform(self):
        m = self.m
        strainrate = 1.0 / m.strainStep * (m.maxStrain / abs(m.maxStrain))
        runTime = m.maxStrain / strainrate / m.timestep
        totalTime = 0
        import StringIO
        f = StringIO.StringIO()
        print >>f, "fix  3  all deform 1 x erate %f remap x units box" % strainrate
        print >>f, "run %d" % (runTime)
        totalTime += runTime
        if m.reverseStrain:
            print >>f, "unfix 3"
            print >>f, "fix  3  all deform 1 x erate %f remap x units box" % (
                -strainrate / (1.0 + m.maxStrain))
            print >>f, "run %d" % (runTime)
            totalTime += runTime
            print >>f, "unfix 3"
            strainrate = 1.0 / m.strainStep * (m.minStrain / abs(m.minStrain))
            runTime = m.minStrain / strainrate / m.timestep
            print >>f, "fix  3  all deform 1 x erate %f remap x units box" % (
                strainrate)
            print >>f, "run %d" % (runTime)
            totalTime += runTime
            print >>f, "unfix 3"
            print >>f, "fix  3  all deform 1 x erate %f remap x units box" % (
                -strainrate / (1.0 + m.minStrain))
            print >>f, "run %d" % (runTime)
            totalTime += runTime
        content = f.getvalue()
        f.close()
        return totalTime, content

    def vDeform(self):
        m = self.m
        # initial length
        lx = m.box[6]

        strainrate = 1.0 / m.strainStep * (m.maxStrain / abs(m.maxStrain))
        runTime = m.maxStrain / strainrate / m.timestep
        totalTime = 0
        import StringIO
        f = StringIO.StringIO()
        print >>f, "fix 3 stayr nve"
        print >>f, "fix force stayr setforce  0 0 0"
        print >>f, "velocity stayr set %f 0.0 0.0 units box" % (
            strainrate * lx)
        print >>f, "run %d" % (runTime)
        totalTime += runTime
        if m.reverseStrain:
            print >>f, "velocity stayr set %f 0.0 0.0 units box" % (
                -strainrate * lx)
            print >>f, "run %d" % (runTime)
            totalTime += runTime
            strainrate = 1.0 / m.strainStep * (m.minStrain / abs(m.minStrain))
            runTime = m.minStrain / strainrate / m.timestep
            print >>f, "velocity stayr set %f 0.0 0.0 units box" % (
                strainrate * lx)
            print >>f, "run %d" % (runTime)
            totalTime += runTime
            print >>f, "velocity stayr set %f 0.0 0.0 units box" % (
                -strainrate * lx)
            print >>f, "run %d" % (runTime)
            totalTime += runTime
        content = f.getvalue()
        f.close()
        return totalTime, content

    def post(self):
        import pandas as pd
        df = pd.read_csv(
            "strain_stress.txt",
            sep=r"[ \t]",
            engine="python",
            skiprows=2,
            header=None)
        import numpy as np
        df = np.array(df)
        strain = df[:, 1] / df[0, 1] - 1
        # convert to GPa
        # metalE/metalV/(siE/siV)=metalE/siE*siV/metalV=1/u.si.E()*u.si.L()**2
        tl.cd('minimize')
        self.m.postMini()
        tl.cd('..')
        pxx = -df[:, 2] * 1e5 * 1e-9 * self.m.zfactor
        from aces.graph import plot
        plot((strain, 'Strain'), (pxx, 'Stress (GPa)'),
             'stress.png', linewidth=2, grid=True)
        dp = pd.DataFrame()
        dp['Strain'] = strain
        dp['Stress_GPa'] = pxx
        dp.to_csv('cal_stress.txt', sep='\t', index=False, float_format="%f")
        from scipy import stats
        slope, intercept, r_value, p_value, std_err = stats.linregress(strain[
                                                                       :50], pxx[:50])
        self.YoungsModulus = slope
        tl.write(slope, 'YoungsModulus.txt')
