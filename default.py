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
	latx=1,laty=1,latz=1,
	ylen=units.metal.L(20),
	thick=units.metal.L(3.35),
	deta=units.metal.L(3),
	T=units.metal.T(300),
	dT=units.metal.T(10),
	metropolis=False,
	write_structure=False,
	timestep=units.metal.t(.5e-3),
	equTime=500000,
	langevin=0,
	nvt=True,jprofile=0,computeTc=1,fourierTc=0,gstart=20000,jcf=0,
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
	begin=10,
	masses='',
	potential='',
	dumpxyz=1,
	dumpv=0,
	hdeta=units.metal.L(3),
	runner='mdTc'
	
	,Cinterval=5
	,Ctime=200000

	,supercell=[1,1,1]
	,supercell3=[1,1,1]
	,useS3=False
	,kpoints=[20,20,1]
	,ekpoints=[3,3,1]
	,ecut=400
	,paw=True
	,gga=True
	,pbe=True
	,mekpoints=[25,25,1]
	,conti=False
	,engine="lammps"
	
	,correlation_supercell=[1,1,1]
	,pho3bte=False
	,corrNVT=False
	,strainStep=5000
	,maxStrain=0.4
	,reverseStrain=False
	,minStrain=-0.4
	,vStrain=True
	,copyN=-1
	,shengcut=-4
	,phofc=False
	,leads='graphene'
	,leadlat=[2,1,1]
	,atomfile=0
	,creatbonds=-1.0
	,phanaonly=False
	,usephana=False
	,dimension=3
	,nseed=10
	,isym=True
	,nkf=[10,10,10]
	,nqf=[10,10,10]
	,ismear=0
	,th=False
)
default['Thi']=default['T']+default['dT']
default['Tlo']=default['T']-default['dT']

	


