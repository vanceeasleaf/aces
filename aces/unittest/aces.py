# encoding: utf-8
# @the host class for a whole project in Automatical Computational Experimental System
# @author Yang Zhou @Fudan University 2015.4.3
import sys,os
class Runner:
	
	# process the consol command and excute the host program
	# from the point of TDD we have to seperate this method from Aces
	def run(aces):
		if sys.argc==1:
			# no parameters so we setup a new experiment
			aces.runSubmit()
		elif sys.arg[1]=='q':
			aces.runQuery()
		elif sys.arg[1]=='clean':
			aces.clean()	
		elif sys.arg[1]=='stop':
			aces.stop()
		else :
			print 'Unkown command!'
			sys.exit() 
			
# there are some strategies to group the apps .
# 1、each app is a group and they have the same cores
# every group uses queue q1.4 's 1 nodes* 4 cores for each app
class NormalStrategy:
	def __init__(self,procs=('q1.4',1,4)):
		self.procs=procs
		self.heaters=[]
		
	def prepare(self,apps):
		self.apps=apps
		for app in apps:
			#mkdir of it
			app.prepare()
			heater=PbsHeater()
			heater.register(app)
			self.heaters.append(heater)
			# prepare pbs files
			heater.prepare(self.procs)
			
# 2、apps are grouped into less groups. each group has the same cores, for instance, 84 cores are shared by 12 apps , and at least 24 cores are accepted as a task,
# so we group them into 3, and each group has 4*7=28 cores, 
# every group uses queue q1.2 's 12 nodes * 2 cores ,per app need 4 cores	
class UniverseStrategy:		
	def __init__(self,procs=('q1.2',12,2,4)):	
		self.procs=procs
		self.heaters=[]
		
	def prepare(self,apps):
		self.apps=apps
		queue,nodes,cores,appCores=self.procs
		appsPerGroup=nodes*cores/appCores
		i=0
		for app in apps:
			#mkdir of it
			app.prepare()
			if i%appsPerGroup==0:
				heater=PbsHeater()
				self.heaters.append(heater)
			heater.register(app)
			heater.prepare(self.procs)
		
# 3、all the apps are in a same group , for example , they run in a single computer.
# using the single computer's 4 cores
class SingleStrategy:
	def __init__(self,procs=(4)):
		self.procs=procs

class VariableStrategy:
	def __init__(self,procs=(84,'q1.2'),target='natom'):
		self.procs=procs			
				
class Aces:
	def __init__(self):
		pass
	
	# the default option of this class
	def setup(self):
		self.strategy=NormalStrategy(procs=(1,4,'q1.4'))
	
	# from the json rebuild the whole structure
	def rebuild(self):
		filename='qloop.txt'
		if not os.path.exists(filename):
			print 'rebuild failed!'
			sys.exit()
		f=open(filename)

		# what stratety, it's cores and queue , the groups and their apps	
		self.strategy.rebuild(f)
		# app parameters
		self.apps.rebuild(f)
	
	
	# stop the apps and clean all the files to the inital state
	def clean():
		self.rebuild()
		pass
		
	def stop():
		s=raw_input("Comfirm to stop all the simulation in this project?[y/n]")
		if s!='y':
			print('exit with no change.')
			sys.exit()
		self.rebuild()
		self.comfirmStop()
		
			
	def comfirmStop():
		pass
		
	# rebuild the apps and excute their query
	def runQuery(self):
		
		#override the original parameter
		self.query()
		self.rebuild()

		#query each app
		self.strategy.query()
			
	# from structure record the json
	def record(self):
		f=open('qloop.txt','w')
		# what stratety, it's cores and queue , the groups and their apps	
		self.strategy.record(f)
		# app parameters
		self.apps.record(f)

	
	def register(self,app):
		self.apps.append(app)
	
	# prepare the apps and excute them
	def runSubmit(self):
		self.submit()
		# assign apps to groups, and preprocess
		self.strategy.prepare(self.apps)
		# excute
		self.strategy.run()
		self.record()
		
	# called before runQuery to setup the conditions, need to be override
	def query(self):
		pass
	
	# called before runSubmit to setup the conditions, need to be override	
	def submit(self):
		pass