# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2017-06-18 21:53:21
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-18 21:51:11

from __future__ import print_function
import os
import sys
import subprocess as sub
import json

printCommand = True


def loadjson(file):
    return json.loads(read(file))


def shell_exec(cmd):
    if printCommand:
        print("[Command]" + cmd)
        sys.stdout.flush()
    c = os.popen(cmd).read()
    return c.strip()


def toString(m, sep=' '):
    return sep.join(map(str, m))


def qdel(a):
    x = "qdel " + toString(a)
    passthru(x)


def passthru(cmd):
    # print os.popen(cmd).read()
    # sys.stdout.flush()
    # sub.call(shlex.split(cmd),stdout=sys.stdout)
    if printCommand:
        print("[Command]" + cmd)
    sys.stdout.flush()
    sub.call(cmd, shell=True, stdout=sys.stdout)


def write(cmd, fileName, mode="w", sep=""):
    file = open(fileName, mode)
    file.write(str(cmd) + sep)
    file.close()


def debug(cmd):
    write(str(cmd), 'debug.txt', 'a', sep="\n")


def sleep(u):
    os.sleep(u)


def exists(path):
    return os.path.exists(path)


def read(fileName):
    file = open(fileName)
    s = file.read()
    file.close()
    return s


def to_txt(columns, data, filename):
    import pandas as pd
    quants = pd.DataFrame(data, columns=columns)
    quants.to_csv(filename, sep='\t', index=None)


def pwrite(fp, s):
    print(s, end="")
    fp.write(s)
    sys.stdout.flush()


def exit(info='Exited by user!'):
    print(info)
    sys.stdout.flush()
    sys.exit()


def mkcd(path):
    mkdir(path)
    cd(path)


def cd(path):
    os.chdir(path)


def mv(src, dest):
    shell_exec("mv %s %s" % (src, dest))


def cp(src, dest):
    shell_exec("cp %s %s -r " % (src, dest))


def mkdir(path):
    path = path.strip()
    path = path.rstrip("\\")
    isExists = os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        return True
    else:
        return False


def ls(path='*'):
    import glob
    return glob.glob(path)


def pwd():
    return os.getcwd()


def dirname(path):
    return os.path.dirname(path).strip()


def basename(path):
    return os.path.basename(path).strip()


def parseyaml(filename):
    try:
        import yaml
    except ImportError:
        print("You need to install python-yaml.")
        exit(1)

    try:
        from yaml import CLoader as Loader
    except ImportError:
        from yaml import Loader
    string = open(filename).read()
    data = yaml.load(string, Loader=Loader)
    return data
