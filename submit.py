#encoding:utf8

import sys,os
from query import clean,stop,query
from toolsub import toolsub
class ACES:
	def __init__(self):
		self.idx=0
		self.projHome=os.path.dirname(sys.argv[0])
		self.projName=os.path.basename(self.projHome);
		self.single=''
	def run(self):
		projHome=self.projHome
		projName=self.projName
		single=self.single
		if(len(sys.argv)==1):
			clean(projHome,projName,single)
			self.subq()
			sys.exit();
		elif(sys.argv[1]=="q"):
			query(projHome,srcHome,universe)
		elif(sys.argv[1]=="clean"):
			clean(projHome,projName,single)
		elif(sys.argv[1]=="stop"):
			stop(projHome,projName,single)
		else:
			print("Unkown command!");
			pass
	def submit(self,opt,cmd):
		projHome=self.projHome
		projName=self.projName
		single=self.single
		species=opt['species']
		units=opt['units']
		method=opt['method'];queue=opt['queue'];nodes=opt['nodes'];procs=opt['procs'];
		toolsub(cmd,self.idx,projHome,projName,species,units,method,queue ,nodes ,procs,universe='',uqueue='',single='',unodes='',uprocs='')
		self.idx+=1
		
	def subq(self):
		pass