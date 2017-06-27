# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   1970-01-01 08:00:00
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-27 22:41:33
import aces.tools as tl
import numpy as np


def meshyaml(file="mesh.yaml"):
    freqs = []
    eigs = []
    qpoints = []
    data = tl.parseyaml(file)
    for phonon in data['phonon']:
        qp = phonon['q-position']
        qpoints.append(qp)
        freq_q = []
        eig_q = []
        for band in phonon['band']:
            frequency = band['frequency']
            freq_q.append(frequency)
            eig = []
            for at in band['eigenvector']:
                eig.append(at)
            eig_q.append(eig)
        eigs.append(eig_q)
        freqs.append(freq_q)
    # shape=(nqpoints,nbranch)
    freqs = np.array(freqs)

    # shape=(nqpoints,nbranch,natom,3,2)
    eigs = np.array(eigs)

    # shape=(nqpoints,nbranch,natom,3) where eigenvector become to complex
    eigs = eigs[..., 0] + 1j * eigs[..., 1]

    # shape=(nqpoints,3)
    qpoints = np.array(qpoints)

    return {'frequency': freqs, "eigs": eigs, "qpoints": qpoints}


def get_eigs(file="mesh.npz"):
    if not tl.exists(file):
        name = file.replace('npz', 'yaml')
        print("generate " + file)
        data = meshyaml(name)
        np.savez(file, **data)
    data = np.load(file)
    freqs = data['frequency']
    eigs = data['eigs']
    qpoints = data['qpoints']
    return qpoints, freqs, eigs


def get_pr(eigs):
    """participation ratio calculation

    $Pr_{k\sigma}= \frac{1}{N \sum_{i\alpha}
    {(\epsilon^{\ialpha}_{k\sigma} \cdot
     \epsilon^{*\ialpha}_{k\sigma})^2}}$

    Arguments:
        eigs {array[nqpoints,nbranch,natom,3]} -- eigenvector
    """

    shape = eigs.shape
    nqpoints, nbranch, natom, n3 = shape[-2]
    eigs = eigs.reshape([nqpoints, nbranch, -1])
    eigs = (eigs * eigs.conjugate())
    eigs = (eigs * eigs).sum(axis=2)

    return 1.0 / eigs / natom
