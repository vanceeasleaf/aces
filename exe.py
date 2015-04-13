#encoding:utf8
from os.path import *
import sys
import os,json,imp
home=dirname(realpath(__file__))
#app home 
projHome=dirname(realpath(sys.argv[0]))
f=open('app.json')
opt=f.read()
opt=json.loads(opt)
species=opt['species']
m= imp.load_source('structure', home+'/materials/'+species+'/structure.py') 
m=m.structure(home,opt)
input= imp.load_source('input', home+'/input.py') 
m.dtime=m.timestep*100;
m.excNum=m.aveRate/m.excRate;
m.corNum=m.aveRate/m.corRate;
m.swapEnergyRate=m.swapEnergy/(m.excRate*m.timestep);
if(m.method=="nvt"):m.xp=0;
m.tcfactor=15;#the unit is J/s/m/K
m.kb=1.0
m.nktv=1.0
input.input(m.units ,m.xp ,m.yp ,m.zp ,m.dumpRate ,m.timestep ,m.method ,m.kb ,m.nktv ,m.masses,m.potential ,m.T ,m.seed ,m.dtime ,m.equTime ,
m.langevin ,m.nvt ,m.aveRate ,m.deta ,m.jprofile  ,m.corRate ,m.computeTc  ,m.fourierTc ,m.tcfactor ,m.gstart ,m.jcf  ,m.nswap ,m.excRate  ,
m.excNum ,m.swapEnergyRate ,m.dumpxyz ,m.dumpv ,m.runTime,m.upP,m.wfix,m.nstat,m.enforceThick,m.thick,m.Thi,m.Tlo,m.hdeta)



	


