constants={
'real':{
'boltz': 0.0019872067,
'hplanck': 95.306976368,
'mvv2e': 48.88821291 * 48.88821291,
'ftm2v': 1.0 / 48.88821291 / 48.88821291,
'mv2d': 1.0 / 0.602214179,
'nktv2p': 68568.415,
'qqr2e': 332.06371,
'qe2f': 23.060549,
'vxmu2f': 1.4393264316e4,
'xxt2kmu': 0.1,
'e_mass': 1.0/1836.1527556560675,
'hhmrr2e': 0.0957018663603261,
'mvh2r': 1.5339009481951,
'angstrom': 1.0,
'femtosecond': 1.0,
'qelectron': 1.0,
'kelvin':1,
'timestep': 1.0,
'ar_mass':39.95
},
'metal':{
'boltz': 8.617343e-5,
'hplanck': 4.135667403e-3,
'mvv2e': 1.0364269e-4,
'ftm2v': 1.0 / 1.0364269e-4,
'mv2d': 1.0 / 0.602214179,
'nktv2p': 1.6021765e6,
'qqr2e': 14.399645,
'qe2f': 1.0,
'vxmu2f': 0.6241509647,
'xxt2kmu': 1.0e-4,
'e_mass': 0.0,    # not yet set
'hhmrr2e': 0.0,
'mvh2r': 0.0,
'angstrom': 1.0,
'femtosecond': 1.0e-3,
'qelectron': 1.0,
'kelvin':1,
'timestep':0.001,
'ar_mass':39.95,
},
"si":{
'boltz': 1.3806504e-23,
'hplanck': 6.62606896e-34,
'mvv2e': 1.0,
'ftm2v': 1.0,
'mv2d': 1.0,
'nktv2p': 1.0,
'qqr2e': 8.9876e9,
'qe2f': 1.0,
'vxmu2f': 1.0,
'xxt2kmu': 1.0,
'e_mass': 0.0,    # not yet set
'hhmrr2e': 0.0,
'mvh2r': 0.0,
'angstrom': 1.0e-10,
'femtosecond': 1.0e-15,
'qelectron': 1.6021765e-19,
'kelvin':1,
'timestep':1.0e-8,
'ar_sigma':3.405e-10,#in si
'ar_epsilon':1.67e-21,
'ar_mass':6.633e-26
},
"cgs":{
'boltz': 1.3806504e-16,
'hplanck': 6.62606896e-27,
'mvv2e': 1.0,
'ftm2v': 1.0,
'mv2d': 1.0,
'nktv2p': 1.0,
'qqr2e': 1.0,
'qe2f': 1.0,
'vxmu2f': 1.0,
'xxt2kmu': 1.0,
'e_mass': 0.0,  # not yet set
'hhmrr2e': 0.0,
'mvh2r': 0.0,
'angstrom': 1.0e-8,
'femtosecond': 1.0e-15,
'qelectron': 4.8032044e-10,
'kelvin':1,
'timestep':1.0e-8,
'ar_mass':6.633e-23
},
"electron":{
'boltz': 3.16681534e-6,
'hplanck': 0.1519829846,
'mvv2e': 1.06657236,
'ftm2v': 0.937582899,
'mv2d': 1.0,
'nktv2p': 2.94210108e13,
'qqr2e': 1.0,
'qe2f': 1.94469051e-10,
'vxmu2f': 3.39893149e1,
'xxt2kmu': 3.13796367e-2,
'e_mass': 0.0,   # not yet set
'hhmrr2e': 0.0,
'mvh2r': 0.0,
'angstrom': 1.88972612,
'femtosecond': 0.0241888428,
'qelectron': 1.0,
'kelvin':1,
'timestep':0.001,
'ar_mass':39.95e-3
},
"lj":{
'boltz': 1.0,# J/K
'hplanck': 0.18292026, # using LJ parameters for argon
'mvv2e': 1.0,
'ftm2v': 1.0,
'mv2d': 1.0,
'nktv2p': 1.0,
'qqr2e': 1.0,
'qe2f': 1.0,
'vxmu2f': 1.0,
'xxt2kmu': 1.0,
'e_mass': 0.0,    # not yet set
'hhmrr2e': 0.0,
'mvh2r': 0.0,
'ar_mass':1.0,
'qelectron': 1.0,
'timestep':0.005
}
}
from math import sqrt;
constants['lj']['angstrom']=constants["si"]['angstrom']/constants["si"]['ar_sigma'];
constants['lj']['femtosecond']=constants["si"]['femtosecond']*sqrt(constants["si"]['ar_epsilon']/constants["si"]['ar_mass']/constants["si"]['ar_sigma']/constants["si"]['ar_sigma'])
constants['lj']['kelvin']=constants["si"]['kelvin']*(constants["si"]['boltz']/constants["si"]['ar_epsilon'])
class Units:
	def __init__(self,targetUnit):
		self.metal=unitBase(targetUnit,'metal')
		self.lj=unitBase(targetUnit,'lj')
		self.real=unitBase(targetUnit,'real')
		self.si=unitBase(targetUnit,'si')
		self.electron=unitBase(targetUnit,'electron')
		self.cgs=unitBase(targetUnit,'cgs')
		# current unit of environment
		self.targetUnit=targetUnit
		units=targetUnit
		si=self.si
		self.tcfactor=1/si.E()*si.t()*si.L()*si.T()
	def __str__(self):
		return self.targetUnit
	def __getattr__( self, name ):
		return constants[self.targetUnit][name]
class unitBase:
	def __init__(self,targetUnit,srcUnits):
		self.targetUnit=targetUnit
		self.srcUnits=srcUnits
	def t(self,val=1.0):
		return  val*self.unitConvert('t')
	def L(self,val=1.0):
		return  val*self.unitConvert('L')	
	def T(self,val=1.0):
		return  val*self.unitConvert('T')
	def M(self,val=1.0):
		return  val*self.unitConvert('M')
	def E(self,val=1.0):
		return  val*self.unitConvert('E')
	def unitConvert(self,type):
		units=self.targetUnit
		srcUnits=self.srcUnits
		if type=='L':return constants[units]['angstrom']/constants[srcUnits]['angstrom']
		if type=="t":return constants[units]['femtosecond']/constants[srcUnits]['femtosecond']
		if type=="T":return constants[units]['kelvin']/constants[srcUnits]['kelvin']
		if type=="M":return constants[units]['ar_mass']/constants[srcUnits]['ar_mass']
		if type=="E":
			a=constants[units]['boltz']*constants[units]['kelvin']
			b=constants[srcUnits]['boltz']*constants[srcUnits]['kelvin']
			return a/b


