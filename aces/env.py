# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2017-06-18 22:01:52
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-18 22:02:30
from os.path import realpath, dirname, basename
from aces.tools import pwd, exists
SRCHOME = dirname(realpath(__file__))


def checkParent(dir, n=0):
    if n == 5:
        raise Exception('error when find sub.py')
    if not exists(dir + '/sub.py'):
        dir = realpath(dir + '/..')
        return checkParent(dir, n + 1)
    else:
        return dir
PROJHOME = checkParent(pwd())
PROJNAME = basename(PROJHOME)
