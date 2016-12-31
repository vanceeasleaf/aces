# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2016-12-25 21:40:33
# @Last Modified by:   YangZhou
# @Last Modified time: 2016-12-25 21:10:07
import sys
from aces.tools import *
print sys.argv
if(len(sys.argv)==1):exit("no arguments");
s=sys.argv[1]
print "Comfirm to kill %s?[y/n]"%s,
q=raw_input(); 
if(q!="y"):exit("exit with no change.")
print "killing ",s
passthru("sudo deljobs.sh %s"%s)