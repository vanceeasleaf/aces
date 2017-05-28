#from aces.env import SRCHOME,PROJHOME,PROJNAME
import json,os.path,time
from aces.runners.minimize import minimize as minimize_input
from importlib import import_module as im
from aces import profile,config
from aces.tools import *
class App:
	def __init__(self):
		"""
		while not os.path.exists('app.json'):
			time.sleep(1)
			print pwd()+'/app.json'
		"""
		if exists('app.json'):
			opt=loadjson('app.json')
		elif exists('../app.json'):
			opt=loadjson('../app.json')
		elif exists('../../app.json'):
			opt=loadjson('../../app.json')
		elif exists('../../../app.json'):
			opt=loadjson('../../../app.json')
		else:exit('app.json lost!')
		species=opt['species']
		s=im('aces.materials.%s'%species)
		m=s.structure(opt)
		self.m=m
		m.home=pwd()
		assert m.home!=''

		Runner=im('aces.runners.%s'%m.runner)
		self.runner=Runner.runner(m)
		
	def minimize(self):	
		if(self.m.copyN==-1):copymini=False
		else: copymini=True
		if copymini:
			while not exists('../%d/minimize/done'%self.m.copyN):
				sleep(30)
			print 'copymini'
			cp('../%d/minimize'%self.m.copyN,'.')
		else:self.creatmini()
				
	def creatmini(self):
		print 'creatmini'
		mkdir('minimize')
		cd('minimize')
		minimize_input(self.m)
		write(time.strftime('%Y-%m-%d %H:%M:%S'),'done')
		cd('..')
	def execute(self):		
		
		self.minimize()
		self.runner.run()

		
	
	def result(self):
		m=self.m
		cd('minimize')
		m.postMini()
		cd('..')
		profile.run(**m.__dict__)

class Apps:
	def __init__():
		obj=[json.loads(json_string) for  json_string in open("qloops.txt")]
if __name__=='__main__':
	App().execute()