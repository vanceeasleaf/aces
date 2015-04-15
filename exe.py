#encoding:utf8
from os.path import *
import sys
import os,json,imp
home=dirname(realpath(__file__))
from aces.Units import Units
from aces.input import input
#app home 
projHome=dirname(realpath(sys.argv[0]))
f=open('app.json')
opt=f.read()
opt=json.loads(opt)
f.close()
species=opt['species']
m= imp.load_source('structure', home+'/materials/'+species+'/structure.py') 
m=m.structure(home,opt)
units=Units(m.units)
m.dtime=m.timestep*100;
m.excNum=m.aveRate/m.excRate;

m.swapEnergyRate=m.swapEnergy/(m.excRate*m.timestep);
if(m.method=="nvt"):m.xp=0;
m.tcfactor=units.tcfactor;#the unit is J/s/m/K
m.kb=units.boltz
m.nktv=units.nktv2p
fixud=m.fixud
input(m.units ,m.xp ,m.yp ,m.zp ,m.dumpRate ,m.timestep ,m.method ,m.kb ,m.nktv ,m.masses,m.potential ,m.T ,m.seed ,m.dtime ,m.equTime ,
m.langevin ,m.nvt ,m.aveRate ,m.deta ,m.jprofile  ,m.corRate ,m.computeTc  ,m.fourierTc ,m.tcfactor ,m.gstart ,m.jcf  ,m.nswap ,m.excRate  ,
m.excNum ,m.swapEnergyRate ,m.dumpxyz ,m.dumpv ,m.runTime,m.upP,m.wfix,m.nstat,m.enforceThick,m.thick,m.Thi,m.Tlo,m.hdeta,fixud)



	


