# encoding: utf-8
# @the host class for a whole project in Automatical Computational Experimental System
# @author Yang Zhou @Fudan University 2015.4.3
import sys,os,json
from aces.tools import exit
#from aces.env import *
class Aces:
	def logo(self):
		print """    _    ____ _____ ____  
   / \  / ___| ____/ ___| 
  / _ \| |   |  _| \___ \ 
 / ___ \ |___| |___ ___) |
/_/   \_\____|_____|____/ 
"""
		print "\nWelcome to Automatical Computational Experiment System(ACES)"
		print "developed by Yang Zhou @Fudan University\n"
		#print "Project Home="+PROJHOME
		#print "Project Name="+PROJNAME
		print ""		
	def __init__(self):
		self.idx=0
		self.logo()
		self.single=''
	def run(self):

		from aces.query import clean,stop,query
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
		from aces.toolsub import toolsub
		origin=dict(single=False
		,species='graphene'
		,units='metal'
		,method='nvt'
		,queue='q1.1'
		,nodes=1
		,procs=4
		,runTime=10000000
		,universe=False
		,uqueue='q1.2'
		,unodes=12
		,uprocs=2
		,runner='mdTc'
		)
		opt=dict(origin,**opt)
		app=dict(opt,**app)
		toolsub(self.idx,app)
		self.idx+=1
		
	# called before runQuery to setup the conditions, need to be override
	def query(self):
		pass
	
	# called before runSubmit to setup the conditions, need to be override	
	def submit(self):
		pass