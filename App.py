from aces.env import SRCHOME,PROJHOME,PROJNAME
import json
from aces.apps.minimize import input as minimize_input
from importlib import import_module as im
from aces import profile,config
from aces.tools import *
class App:
	def __init__(self):
		f=open('app.json')
		opt=f.read()
		opt=json.loads(opt)
		f.close()
		species=opt['species']
		s=im('aces.materials.%s.structure'%species)
		m=s.structure(opt)
		self.m=m
		
	def minimize(self):	
		mkdir('minimize')
		cd('minimize')
		minimize_input(self.m)
		shell_exec(config.mpirun+" %s "self.m.cores+config.lammps+" <input >log.out")
		cd('..')
		
	def execute(self):
		self.minimize()
		self.post()
		m=self.m
		Runner=im('aces.runners.%s'%m.runner)
		runner=Runner.runner(m)
		runner.run()

		
	def post(self):
		self.m.postMini()
	
	def result(self):
		m=self.m
		self.post()
		profile.run(**m.__dict__)

class Apps:
	def __init__():
		obj=[json.loads(json_string) for  json_string in open("qloops.txt")]
if __name__=='__main__':
	App().execute()