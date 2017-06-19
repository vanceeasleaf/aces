# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2017-06-19 19:10:08
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-19 19:10:14
import sys
from aces.tools import shell_exec


class Runner:

    def __init__(self, m):
        self.m = m

    def generate(self):
        pass

    def run0(self):
        __console__ = sys.stdout

        f = open('input', 'w', 0)
        sys.stdout = f
        self.generate()
        sys.stdout = __console__
        f.close()
        cmd = self.runcmd()
        shell_exec(cmd)

    def runcmd(self):
        return ""

    def post(self):
        pass
