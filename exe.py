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
input(m)



	


