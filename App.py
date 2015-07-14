from aces.env import SRCHOME,PROJHOME,PROJNAME
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
		f=open('app.json')
		opt=f.read()
		opt=json.loads(opt)
		f.close()
		species=opt['species']
		s=im('aces.materials.%s'%species)
		m=s.structure(opt)
		self.m=m
		Runner=im('aces.runners.%s'%m.runner)
		self.runner=Runner.runner(m)
		
	def minimize(self):	
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