#encoding:utf8
from aces.env import *
import os,json,imp

from aces.input import input

f=open('app.json')
opt=f.read()
opt=json.loads(opt)
f.close()
species=opt['species']
m= imp.load_source('structure', SRCHOME+'/materials/'+species+'/structure.py') 
m=m.structure(SRCHOME,opt)
input(m)



	


