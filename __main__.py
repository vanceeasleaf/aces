# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2016-09-05 19:33:20
# @Last Modified by:   YangZhou
# @Last Modified time: 2016-09-05 19:31:44

from aces.App import App
from aces.tools import *
import aces.tools
import sys
import json
def run():
	aces.tools.printCommand=False
	r=App().runner
	if len(sys.argv)>1:
		a=sys.argv[1]
		if not hasattr(r,a):
			print a,'方法不存在'
		else: getattr(r,a)()
if exists("sub.py"):
	obj=[json.loads(json_string) for  json_string in open("qloops.txt")]
	for i in range(len(obj)):
		l=str(i)
		cd(l)
		run()
		cd('..')
else:
	run();
