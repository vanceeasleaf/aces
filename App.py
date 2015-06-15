from aces.env import SRCHOME,PROJHOME,PROJNAME
import json
from aces.minimize.input import input as minimize_input
from aces.input import input as exe_input
from aces.input import postMini
from importlib import import_module as im
from aces.Units import Units
from aces import profile
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
		m=self.m
		minimize_input(m)
	def execute(self):
		m=self.m
		self.post()
		exe_input(m)
	def post(self):
		m=self.m
		units=Units(m.units)
		m.kb=units.boltz
		m.nktv=units.nktv2p
		if(m.method=="nvt"):m.xp=0;
		lx,ly,lz,m.zfactor,m.S,xlo,xhi,ylo,yhi,zlo,zhi=postMini(m.xp,m.yp,m.zp,m.enforceThick,m.thick)
		m.dtime=m.timestep*100;
		m.tcfactor=units.tcfactor;
		m.excNum=m.aveRate/m.excRate;
		m.swapEnergyRate=m.swapEnergy/(m.excRate*m.timestep);
		m.units=units
		m.box=(xlo,xhi,ylo,yhi,zlo,zhi,lx,ly,lz)		
	def result(self):
		m=self.m
		self.post()
		profile.run(**m.__dict__)

class Apps:
	def __init__():
		obj=[json.loads(json_string) for  json_string in open("qloops.txt")]
		