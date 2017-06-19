# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2017-06-19 13:11:33
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-19 13:11:56
from aces._io.lineManager import lineManager


class strain:

    def __init__(self, dumpname):
        self.lm = lineManager(dumpname)
        self.info()

    def info(self):
        lm = self.lm
        self.natom = int(lm.getLine(3).split()[0])
        t1 = int(lm.getLine(1).split()[0])
        self.line_interval = 9 + self.natom
        if lm.nline < self.line_interval:
            self.interval = 1
        else:
            t2 = int(lm.getLine(1 + self.line_interval).split()[0])
            self.interval = t2 - t1
        self.totalStep = lm.nline / self.line_interval
        if self.totalStep % 2 == 1:
            self.totalStep -= 1
        print("Atom Number=", self.natom)
        print("Total step=", self.totalStep)
        print("interval=", self.interval)

    def grep(self, index):
        s = ""
        for i in range(self.line_interval):
            s += self.lm.getLine(i + index * self.line_interval) + "\n"
        return s
