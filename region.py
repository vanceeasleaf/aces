# encoding: utf-8
import json
from vapory import *
from random import choice
from ase.io import read
from aces.input import getboxrange
import numpy as np
def drawRegions():
	xlo,xhi,ylo,yhi,zlo,zhi=getboxrange()
	f=open('regions.txt')
	regions=json.loads(f.read())
	f.close()
	objs=[]
	a=[0.99]
	colors=[
	[0,0,1],[0,1,0],[1,0,0],[1,1,0],[1,0,1],[0,1,1],[1,1,1],[2,1,1],[1,2,1],[1,1,2]
	]
	wall = Plane([0, 0, 1], 100,Texture('T_Ruby_Glass'))
	ground = Plane( [0, 1, 0], 0,
                Texture( Pigment( 'color rgb', [1, 1, 1]),
                         Finish( 'phong', 0.8,
                                 'reflection',0.1,
                                 'metallic', 0.9,'ior',1.5,'diffuse', .5)))
	light = LightSource([(xhi+xlo)/2, (yhi+ylo)/2, (zhi+zlo)/2], 'White')
	lo=[xlo,ylo,zlo]
	hi=[xhi,yhi,zhi]
	k=0
	for region in regions:
		if region['type']=='box':
			x=region['dim'][0]
			y=region['dim'][1]
			for i in range(3):
				if x[i]=='INF':x[i]=lo[i]
				if y[i]=='INF':y[i]=hi[i]
			obj = Box(x,y, Pigment('color', 'rgbf',colors[k]+a),Finish('phong', 0,
                                                  'reflection', 0.0,'metallic', .1 ),Interior('ior',1))
			objs.append(obj)
		k+=1	
	scene = Scene( Camera("location", [100, 100, -100], "look_at", [(xhi+xlo)/2, (yhi+ylo)/2, (zhi+zlo)/2],'direction',[0,0,1.5]),objects = [ wall,ground,light]+objs,included=["glass.inc","colors.inc"] )
	scene.render('regions.png',width=800,height=600 ,antialiasing=.01)