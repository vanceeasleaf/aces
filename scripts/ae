# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2016-09-05 19:33:20
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-05-25 09:32:10


from aces.tools import *

import aces.tools
import sys
import json
def run():
	from aces.App import App
	r=App().runner

	if len(sys.argv)>1:
		a=sys.argv[1]
		if not hasattr(r,a):
			print a,'method does not exist'
		else: getattr(r,a)()
def exe():
	if len(sys.argv)>1:
		a=sys.argv[1]
		if(a.find("git")>0):
			mkdir('.project')
			cd('.project')
			if(not exists('.git')):
				shell_exec("git init")
				write("_gsdata_",".gitignore")
			d=ls()
			for x in d:
				if(x.find(".git")>=0):continue
				if(x.find("_gsdata_")>=0):continue
				shell_exec("rm %s -r"%x)
			cd('..')
			passthru("find .  -path ./.project  -prune -o -name '*.py' -print| cpio -pdm .project 2>/dev/null ;cd .project;git add .;git commit -am \"%s\""%(sys.argv[2]))
			return
	from aces.App import App
	from aces import config
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
			elif len(sys.argv)==4 and sys.argv[3]=="-pbs":
				from aces.jobManager import pbs
				job=pbs(queue='q1.1',nodes=4,procs=12,disp="script",path=pwd(),content=config.mpirun+ " 48 "+config.python+' '.join(sys.argv).replace('-pbs','')+' >aces.out')
				cd(job.path)
				job.writepbs()
				job.submit()
			else:getattr(r,a)()
			return
	run();
exe()