#encoding:utf8
import sys,os
sys.path.append(os.path.realpath(os.path.dirname(__file__)))
from Units import Units
units=Units('metal')
default={
	'units':"metal",
	'method':"nvt",
	'enforceThick':1,
	'xp':1,
	'yp':1,#yPeriodic
	'zp':1,
	'useMini':1,
	'dumpxyz':1,
	'upP':3,

	'usinglat':1,
	'latx':4,'laty':4,'latz':4,
	'ylen':units.metal.L(20),
	'thick':units.metal.L(3.35),
	'deta':units.metal.L(3),
	'T':units.metal.T(300),#K for all units

	'metropolis':0,
	'write_structure':0,
	'timestep':units.metal.t(.5e-3),
	'equTime':100000,
	'langevin':0,
	'nvt':1,'jprofile':0,'computeTc':1,'fourierTc':0,'gstart':20000,'jcf':0,
	'dumpRate':100000,
	'aveRate':100000,
	'excRate':500,
	'corRate':2,
	'runTime':10000000,
	'seed':458127641,
	#***********for nvt***********

	'nstat':1,#heat bath width,unit:deta
	'nswap':1,
	'swapEnergy':5e-4,
	'wfix':3,#fix width
	'Nbins':500,
	'bond':units.metal.L(1.42),
	'begin':1,
	'masses':'',
	'potential':'',
	'dumpxyz':1,
	'dumpv':0,
	'hdeta':units.metal.L(3)
}
default['Thi']=default['T']+units.metal.T(10)
default['Tlo']=default['T']-units.metal.T(10)

	


