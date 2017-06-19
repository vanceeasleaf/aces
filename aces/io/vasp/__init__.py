# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2017-06-01 21:49:49
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-19 13:05:06
import numpy as np
import aces.tools as tl
from aces import config


def writevasp(atoms, file='POSCAR'):
    f = open(file, 'w')
    s = np.array(atoms.get_chemical_symbols())
    ss = atoms.get_scaled_positions()
    print >>f, 'ACES POSCAR'
    print >>f, '1.0'
    for x in atoms.cell:
        print >>f, tl.toString(x)
    # ele=np.unique(s)
    ele = []
    for a in s:
        if a in ele:
            continue
        ele.append(a)
    print >>f, tl.toString(ele)
    a = []
    # len(s) = natom
    p = np.arange(len(s))
    for e in ele:
        a.append(p[s == e])
    # p= 0 1 2 3 4 5
    # s= C N N C N C
    # a= [[0,3,5],[1,2,4]]
    ns = [len(x) for x in a]
    # ns =[3,3]
    print >>f, tl.toString(ns)
    print >>f, 'Direct'
    v = []
    for x in a:
        for u in x:
            v.append(u)
            print >>f, tl.toString(ss[u])
    # v= [0,3,5,1,2,4]
    f.close()
    x = np.array(v, dtype=np.int).argsort()
    np.savetxt('POSCARswap', x)


def writePOTCAR(options, elements):
    dir = 'pot'  # LDA
    # paw：PAW-LDA
    # paw_gga：PAW-GGA-PW91
    # paw_pbe：PAW-GGA-PBE
    # pot：USPP-LDA
    # pot_GGA：USPP-GGA
    if not options.paw:
        if options.gga:
            dir = 'pot_GGA'
        else:
            dir = 'pot'
    else:
        if not options.gga:
            dir = 'paw'
        else:
            if options.pbe:
                dir = 'paw_pbe'
            else:
                dir = 'paw_gga'
    tl.passthru('cat "" >POTCAR')
    for ele in elements:
        file = config.vasppot + "/%s/%s/POTCAR" % (dir, ele)
        z = False
        if not tl.exists(file):
            file += '.Z'
            z = True
        assert tl.exists(file)
        if z:
            tl.passthru('zcat %s >> POTCAR' % file)
        else:
            tl.passthru('cat %s >> POTCAR' % file)
    # s=''.join([tools.read(config.vasppot\
    # +"/%s/%s/POTCAR.Z"%(dir,ele)) for ele in self.elements])
    # tools.write(s,'POTCAR')


def parseVasprun(vasprun, tag="forces"):
    collection = []
    for event, element in vasprun:
        if element.attrib['name'] == tag:
            for v in element.xpath('./v'):
                collection.append([float(x) for x in v.text.split()])
    collection = np.array(collection)
    return collection


def writeKPOINTS(kpoints):
    s = """A
    0
    Monkhorst-Pack
    %s
    0  0  0
    """ % ' '.join(map(str, kpoints))
    s = tl.headtrim(s)
    tl.write(s, 'KPOINTS')
