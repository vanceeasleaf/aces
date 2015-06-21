# encoding: utf-8
# @the host class for a whole project in Automatical Computational Experimental System
# @author Yang Zhou @Fudan University 2015.4.3
import sys,os,json
from aces.query import clean,stop,query
from aces.toolsub import toolsub
from aces.tools import exit
from aces.env import *
class Aces:
	def __init__(self):
		self.idx=0
		print """	_    ____ _____ ____  
   / \  / ___| ____/ ___| 
  / _ \| |   |  _| \___ \ 
 / ___ \ |___| |___ ___) |
/_/   \_\____|_____|____/ 
"""
		print "\nWelcome to Automatical Computational Experiment System(ACES)"
		print "developed by Yang Zhou @Fudan University\n"
		print "Project Home="+PROJHOME
		print "Project Name="+PROJNAME
		print ""
		self.single=''
	def run(self):

		single=self.single

		universe=''
		if(len(sys.argv)==1):
			clean(single)
			self.submit()
			sys.exit();
		elif(sys.argv[1]=="q"):
			query(universe)
		elif(sys.argv[1]=="clean"):
			clean(single)
		elif(sys.argv[1]=="stop"):
			stop(single)
		else:
			exit("Unkown command!");
	def commit(self,opt,app):

		single=self.single
		species=opt['species']
		units=opt['units']
		method=opt['method'];
		queue=opt['queue'];
		nodes=opt['nodes'];
		procs=opt['procs'];
		runTime=opt['runTime']
		bte=False
		if opt.has_key('bte'):bte=opt['bte']
		cmd=''
		for key in app:
			val=app[key]
			cmd+='$%s=%s;'%(key,val)
		jj=json.dumps([cmd,app])
		toolsub(cmd,self.idx,species,units,method,queue ,nodes ,procs ,runTime,jj,universe='',uqueue='',single='',unodes='',uprocs='',bte=bte)
		self.idx+=1
		
	# called before runQuery to setup the conditions, need to be override
	def query(self):
		pass
	
	# called before runSubmit to setup the conditions, need to be override	
	def submit(self):
		pass