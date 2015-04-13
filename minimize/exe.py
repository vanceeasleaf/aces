#encoding:utf8
from os.path import *
import sys
import os,json,imp
from aces.minimize.input import input
home=abspath(dirname(realpath(__file__))+'/..');
#app home 
projHome=abspath(dirname(realpath(sys.argv[0]))+'/..')
f=open('../app.json')
opt=f.read()
opt=json.loads(opt)
species=opt['species']
m= imp.load_source('structure', home+'/materials/'+species+'/structure.py') 
m=m.structure(home,opt)
units=opt['units']
input(units,m.structure,m.potential,m.timestep,m.masses,m.dumpRate,m.write_structure,m.metropolis,m.useMini,m.dump)