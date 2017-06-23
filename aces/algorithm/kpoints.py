# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2017-06-23 17:31:09
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-23 17:49:09
"""kpoints related functions"""
import numpy as np


def filter_along_direction(ks, t, eps=0.1 * np.pi, unit="rad"):
    # find the k points that are in the t direction
    if t < 0:
        t = -t
    if unit == "deg":
        t = t * np.pi / 180.0
        eps = eps * np.pi / 180.0
    u = np.arctan2(ks[:, 1], ks[:, 0])
    b = u - t
    b[b > np.pi] -= 2.0 * np.pi
    b[b < -np.pi] += 2.0 * np.pi
    filter1 = np.abs(b) < eps
    return filter1
