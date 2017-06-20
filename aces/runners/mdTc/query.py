# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2017-06-02 18:39:05
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-20 14:57:30

from __future__ import print_function
import os
import aces.config as config
from .inequality import inequality

import aces.tools as tl
from ase.io import read
from aces.env import SRCHOME


def getRatio(path):
    if(not os.path.exists(path)):
        return 0.0
    fp = open(path, "r")
    fp.next()
    natom = int(fp.next().split()[0])
    ntype = int(fp.next().split()[0])
    if ntype == 1:
        return 0.0
    n = 0
    label = ""
    while(label != "Atoms" and n < 20):
        label = fp.next().strip()
        n += 1

    fp.next()
    a = [0.0] * ntype
    for line in fp:
        type = int(line.split()[1])
        a[type - 1] += 1
    return float(a[1]) / natom


def getQueryInfo(workPath, pid, runTime):
    pid = filter(str.isdigit, str(pid))
    lastline = tl.shell_exec("tail -3 %s/log.out" % workPath)
    qstat = tl.shell_exec("qstat %s 2>&1|tail -1 " % pid)

    step = lastline.split()[0]
    if step.isdigit():
        percent = "%.1f%%" % (float(step) / runTime * 100)
    else:
        percent = "0"
    if(qstat.find("Unknown Job Id") >= 0):  # 该任务已从任务队列中去除*/
        time = "complete"
        if(lastline.find("builds") >= 0):
            status = "C"
            percent = "100%"
        else:  # 异常退出*/
            status = "E"

    else:  # 正在运行或等待R&Q&C*/
        time, status, queue = qstat.split()[3:6]
        info = tl.shell_exec("qstat -f %s 2>&1|grep nodes" % pid)
        info = info.split()[2]

    return (percent, status)


def kappa():
    kappaline = tl.shell_exec("tail -1 result.txt 2>err;")
    kappa = kappaline.split('=')
    if len(kappa) > 1:
        kappa = kappa[1]
        return kappa
    return 0.0


def tEnerty():
    """total energy

    [description]

    Returns:
            [type] -- [description]
    """
    totalEline = tl.shell_exec("cd minimize;tail -22 log.out| head -1;")
    a = totalEline.split()
    if(len(a) > 1):
        totalE = [1]
        return float(totalE)
    return 0.0


def nAtom():
    # atom number */
    Natomline = tl.shell_exec("cd minimize;grep atoms log.out ;")
    Natom = Natomline.split()[0]
    if(Natom.isdigit() and Natom > 0):
        return Natom
    return 0


def tDisorder():

    # disorder degree*/
    now = tl.pwd()
    tl.cd('minimize')
    tl.mkcd('disorder')

    disorderLine = tl.shell_exec(
        "cp %s" %
        SRCHOME +
        "/in.disorder .;" +
        config.lammps +
        " <in.disorder 2>err 1>log;tail -1 disorder.txt  2>err;")
    k = disorderLine.split()[1:3]
    if len(k) == 1:
        k.append("")
    disorder, rd = k
    tl.cd(now)
    return (disorder, rd)


def drawStructure():
    atoms = read('minimize/range', format='lammps')
    if len(atoms) < 200:
        atoms.write('minimize.png')


def ineq(m):
    now = tl.pwd()
    species = m.species
    if not (species in ["CN-small"]):
        return 0.0
    tl.cd('minimize')
    tl.mkdir('nonequ')
    tl.cd('nonequ')

    ie = inequality()
    nonequ5 = ie.run()
    tl.cd(now)
    return nonequ5
