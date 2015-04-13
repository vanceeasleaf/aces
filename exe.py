#encoding:utf8
from os.path import *
import sys
import os,json,imp
home=abspath(dirname(realpath(__file__))+'/..');
#app home 
projHome=abspath(dirname(realpath(sys.argv[0]))+'/..')
f=open('../app.json')
opt=f.read()
opt=json.loads(opt)
species=opt['species']
m= imp.load_source('structure', home+'/materials/'+species+'/structure.py') 
m=m.structure(home,opt)
input= imp.load_source('input', home+'/input.py') 
units=opt['units']
method=opt['method']
timestep=opt['timestep']
dtime=timestep*100;
excNum=aveRate/excRate;
corNum=aveRate/corRate;
swapEnergyRate=swapEnergy/(excRate*timestep);
xp=1
if(method=="nvt")xp=0;
tcfactor=15;#the unit is J/s/m/K
input(units ,xp ,yp ,zp ,dumpRate ,timestep ,method ,kb ,nktv ,masses,potential ,T ,seed ,dtime ,equTime ,
langevin ,nvt ,aveRate ,deta ,jprofile ,dumpRate ,corRate ,computeTc  ,fourierTc ,tcfactor ,gstart ,jcf  ,nswap ,excRate  ,
excNum ,swapEnergyRate ,dumpxyz ,dumpv ,runTime,upP,wfix,nstat,enforceThick,thick,Thi,Tlo,hdeta)



	


