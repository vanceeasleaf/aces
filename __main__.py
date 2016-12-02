# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2016-09-05 19:33:20
# @Last Modified by:   YangZhou
# @Last Modified time: 2016-12-02 12:40:58

from aces.App import App
from aces.tools import *
from aces import config
import aces.tools
import sys
import json
def run():
	r=App().runner
	if len(sys.argv)>1:
		a=sys.argv[1]
		if not hasattr(r,a):
			print a,'method does not exist'
		else: getattr(r,a)()
def exe():
	if exists("sub.py"):
		aces.tools.printCommand=False
		obj=[json.loads(json_string) for  json_string in open("qloops.txt")]
		for i in range(len(obj)):
			l=str(i)
			cd(l)
			run()
			cd('..')
		aces.tools.printCommand=True
		return
	if len(sys.argv)>2:
		if sys.argv[2]=="-pbs":
			from aces.jobManager import pbs
			m=App().m
			job=pbs(queue=m.queue,nodes=m.nodes,procs=m.procs,disp=m.pbsname,path=pwd(),content=config.python+' '.join(sys.argv).replace('-pbs','')+' >aces.out')
			cd(job.path)
			job.writepbs()
			job.submit()
			return
		if sys.argv[2]=='-c':
			a=sys.argv[1]
			import aces.script as r
			if not hasattr(r,a):
				print a,'method does not exist'
			else: getattr(r,a)()
			return
	run();
exe()