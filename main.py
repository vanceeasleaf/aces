# encoding: utf-8
# @the host class for a whole project in Automatical Computational Experimental System
# @author Yang Zhou @Fudan University 2015.4.3
import sys,os,json
from aces.query import clean,stop,query
from aces.toolsub import toolsub
class Aces:
	def __init__(self):
		self.idx=0
		self.projHome=os.path.dirname(os.path.realpath(sys.argv[0]))
		self.projName=os.path.basename(self.projHome);
		print "\nWelcome to Automatical Computational Experiment System(ACES)"
		print "developed by Yang Zhou @Fudan University\n"
		print "Project Home="+self.projHome
		print "Project Name="+self.projName
		print ""
		self.single=''
	def run(self):
		projHome=self.projHome
		projName=self.projName
		single=self.single
		srcHome=os.path.dirname(os.path.realpath(__file__))
		universe=''
		if(len(sys.argv)==1):
			clean(projHome,projName,single)
			self.submit()
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
	def commit(self,opt,app):
		projHome=self.projHome
		projName=self.projName
		single=self.single
		species=opt['species']
		units=opt['units']
		method=opt['method'];queue=opt['queue'];nodes=opt['nodes'];procs=opt['procs'];runTime=opt['runTime']
		cmd=''
		for key in app:
			val=app[key]
			cmd+='$%s=%s;'%(key,val)
		jj=json.dumps([cmd,app])
		toolsub(cmd,self.idx,projHome,projName,species,units,method,queue ,nodes ,procs ,runTime,jj,universe='',uqueue='',single='',unodes='',uprocs='')
		self.idx+=1
		
	# called before runQuery to setup the conditions, need to be override
	def query(self):
		pass
	
	# called before runSubmit to setup the conditions, need to be override	
	def submit(self):
		pass