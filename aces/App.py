# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2017-06-18 21:58:57
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-18 21:57:34
#from aces.env import SRCHOME,PROJHOME,PROJNAME
import time
from aces.runners.minimize import minimize as minimize_input
from importlib import import_module as im
import aces.tools as tl


class App:

    def __init__(self):
        """
                while not os.path.exists('app.json'):
                        time.sleep(1)
                        print pwd()+'/app.json'
                """
        if tl.exists('app.json'):
            opt = tl.loadjson('app.json')
        elif tl.exists('../app.json'):
            opt = tl.loadjson('../app.json')
        elif tl.exists('../../app.json'):
            opt = tl.loadjson('../../app.json')
        elif tl.exists('../../../app.json'):
            opt = tl.loadjson('../../../app.json')
        else:
            tl.exit('app.json lost!')
        species = opt['species']
        s = im('aces.materials.%s' % species)
        m = s.structure(opt)
        self.m = m
        m.home = tl.pwd()
        assert m.home != ''

        Runner = im('aces.runners.%s' % m.runner)
        self.runner = Runner.runner(m)

    def minimize(self):
        if (self.m.copyN == -1):
            copymini = False
        else:
            copymini = True
        if copymini:
            while not tl.exists('../%d/minimize/done' % self.m.copyN):
                tl.sleep(30)
            print('copymini')
            tl.cp('../%d/minimize' % self.m.copyN, '.')
        else:
            self.creatmini()

    def creatmini(self):
        print('creatmini')
        tl.mkdir('minimize')
        tl.cd('minimize')
        minimize_input(self.m)
        tl.write(time.strftime('%Y-%m-%d %H:%M:%S'), 'done')
        tl.cd('..')

    def execute(self):

        self.minimize()
        self.runner.run0()


class Apps:

    def __init__():
        pass
