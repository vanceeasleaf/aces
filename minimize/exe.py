#encoding:utf8
from aces.env import SRCHOME
import os,json,imp
from aces.minimize.input import input

f=open('../app.json')
opt=f.read()
opt=json.loads(opt)
species=opt['species']
m= imp.load_source('structure', SRCHOME+'/materials/'+species+'/structure.py') 
m=m.structure(SRCHOME,opt)
units=opt['units']
input(units,m)