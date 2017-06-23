# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   1970-01-01 08:00:00
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-23 20:23:50
import os
import numpy as np


def get_qpoints_full(dir="."):
    ###########################################
    #
    # understanding format of BTE.qpoints_full
    # qpoints_full.shape = (nqpointsf,5)
    # qpoints_full[:,2:] = all qpoints
    # qpoints_full[:,1]=idx
    # qpoints[idx-1]=qpoints_full
    # tau[idx-1]=tau_full
    #
    ###########################################
    file = os.path.join(dir, "BTE.qpoints_full")
    qpoints_full = np.loadtxt(file)

    # the range of qpoints are [0,1],however we want firt BZ [-.5,0.5]
    fil = qpoints_full > 0.5

    # don't change idx columns
    fil[:, :2] = False
    qpoints_full[fil] -= 1.0
    idx = (qpoints_full[:, 1] - 1).astype('int')
    return qpoints_full[:, 2:], idx


def get_qpoints(dir="."):
    ###########################################
    #
    # understanding format of BTE.qpoints
    # qpoints_full.shape = (nqpoints,6)
    # qpoints_full[:,3:] = all qpoints
    # qpoints_full[:,1]=idx
    # qpoints_full[idx-1]=qpoints
    #
    ###########################################
    file = os.path.join(dir, "BTE.qpoints")
    qpoints = np.loadtxt(file)
    return qpoints[:, 3:]


def get_omega(dir="."):
    ###########################################
    #
    # understanding format of BTE.omega
    # omega.shape = (nqpoints,3natom)
    # for q in qpoints:
    #   omega.append([omega of all branchs])
    #
    ###########################################
    file = os.path.join(dir, "BTE.omega")
    omega = np.loadtxt(file) / (2.0 * np.pi)
    return omega


def get_w_final(dir=".", sub="T300K"):
    ###########################################
    #
    # understanding format of BTE.w_final
    # w_final.shape = (N,2), N=3natom*nqpoints
    # for branch in range(3natom)
    #   for q in qpoints:
    #       w_final.append([omega_i,w_i])
    #
    # what we want is (nqpoints,3natom) like omega
    #
    ###########################################
    file = os.path.join(dir, sub, "BTE.w_final")
    w_final = np.loadtxt(file)[:, 1]

    w = np.abs(w_final)
    q = get_qpoints(dir)

    # reshape w to like omega
    n = len(q)

    # (1,N)->(3natom,nqpoints)
    w = w.T.reshape([-1, n])

    # (3natom,nqpoints)
    w = np.einsum('jk->kj', w)
    w.flags.writeable = True

    omega = get_omega(dir)
    # (N,1)
    w[omega < omega.flatten().max() * 0.005] = float('nan')

    return w


def get_tau(dir="."):
    w = get_w_final(dir)
    tau = 1.0 / w + 1e-6
    return tau


def get_v(dir="."):
    ###########################################
    #
    # understanding format of BTE.v
    # v.shape = (N,3), N=3natom*nqpoints
    # for branch in range(3natom)
    #   for q in qpoints:
    #       v.append([v0,v1,v2])
    #
    # what we want is (nqpoints,3natom,3)
    #
    ###########################################
    file = os.path.join(dir, "BTE.v")
    v = np.loadtxt(file)
    q = get_qpoints(dir)
    n = len(q)

    # shape=[3,nbranch,nqpoints]
    v = v.T.reshape([3, -1, n])

    # shape=[nqpoints,nbranch,3]
    v = np.einsum('ijk->kji', v)
    return v


def get_gruneisen(dir="."):
    ###########################################
    #
    # understanding format of BTE.gruneisen
    # gruneisen.shape = (nqpoints,3natom)
    # for q in qpoints:
    #   gruneisen.append([gruneisen of all branchs])
    #
    ###########################################
    file = os.path.join(dir, "BTE.gruneisen")
    gruneisen = np.loadtxt(file)
    return gruneisen
