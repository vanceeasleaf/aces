#encoding:utf8
import sys,os
from aces.Units import Units
units=Units('metal')
default=dict(
	units="metal",
	method="nvt",
	enforceThick=True,
	xp=1,
	yp=1,
	zp=1,
	useMini=True,
	upP=3,
	fixud=0,
	usinglat=True,
	latx=4,laty=4,latz=4,
	ylen=units.metal.L(20),
	thick=units.metal.L(3.35),
	deta=units.metal.L(3),
	T=units.metal.T(300),

	metropolis=False,
	write_structure=False,
	timestep=units.metal.t(.5e-3),
	equTime=500000,
	langevin=0,
	nvt=1,jprofile=0,computeTc=1,fourierTc=0,gstart=20000,jcf=0,
	dumpRate=100000,
	aveRate=100000,
	excRate=500,
	corRate=2,
	runTime=10000000,
	seed=458127641,
	#***********for nvt***********

	nstat=1,#heat bath width,unit:deta
	nswap=1,
	swapEnergy=5e-4,
	wfix=3,#fix width
	Nbins=500,
	bond=units.metal.L(1.42),
	begin=1,
	masses='',
	potential='',
	dumpxyz=1,
	dumpv=0,
	hdeta=units.metal.L(3),
	runner='mdTc'
	
	,Cinterval=5
	,Ctime=200000
	
	,supercell=[1,1,1]
	,kpoints=[20,20,1]
	
	,conti=False
	,engine="lammps"
	
	,correlation_supercell=[1,1,1]
	,pho3bte=False
	,corrNVT=False
	,copymini=False
	,strainStep=5000
	,maxStrain=0.4
	,reverseStrain=False
	,minStrain=-0.4
	,vStrain=False
)
default['Thi']=default['T']+units.metal.T(10)
default['Tlo']=default['T']-units.metal.T(10)

	


