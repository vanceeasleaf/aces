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
	a=[1.0]
	colors=np.array([
	[0,0,1],[0,1,0],[1,0,0],[1,1,0],[1,0,1],[0,1,1],[.9,.5,.9],[2,1,1],[1,2,1],[1,1,2]
	])

	
	wall = Plane([0, 1, 0], -(yhi-ylo)/2,Texture( Pigment( 'color rgb', [1, 1, 1]),
                         Finish( 'phong', 0.8,
                                 'reflection',0.0,
                                 'metallic', 0.1,'ior',1.5,'diffuse', .5)))
	"""
	ground = Plane( [0, 0, -1], -(zhi-zlo)/2,
                Texture( Pigment( 'color rgb', [1, 1, 1]),
                         Finish( 'phong', 0.8,
                                 'reflection',0.7,'ambient',0.5,
                                 'metallic', 0.8,'ior',1.5,'diffuse', .9)))
      """
	ground = Plane( [0, 0, -1], -(zhi+zlo)/2,Texture(Pigment("""    gradient y
      color_map {
         [0,    0.25 color Gray      color Gray]
         [0.25, 0.50 color DimGray   color LightGray]
         [0.50, 0.75 color LightGray color Gray]
         [0.75, 1    color Gray      color Gray]
      }
      scale <1, 8, 1> turbulence 5""")))
	light = LightSource([0, -50,-(zhi-zlo)/2-200], 'White shadowless')
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
			obj = Box(x,y, Pigment('color', 'rgbf',list(1*colors[k])+a),Finish('phong', 0.5,'ambient 0.7',
                                                  'reflection', 0.3,'metallic', 0.2 ),Interior('ior',1.2))
			objs.append(obj)
		k+=1
	atoms=read('minimize/range',format='lammps')
	balls=[]
	for pos in atoms.positions:
		ball=Sphere(np.array(pos)+np.array([xlo,ylo,zlo]),0.7,Pigment('color White', ),Finish('phong', 1,
                                                  'reflection', 0.1,'metallic', .1 ),Interior('ior',1.2)
                                                  )
		balls.append(ball)
	object=Union().add_args(objs+balls+['translate %f*x'%(-(xhi+xlo)/2),'translate %f*y'%(-(yhi+ylo)/2),'translate %f*z'%(-(zhi+zlo)/2)])
	scene = Scene( Camera('orthographic',"location", [0, 0, -2], "look_at", [0, 0,0],'direction',[0,0,1],'sky',[0,0,-1],'scale 70'),objects = [ ground,light,object],included=["glass.inc","colors.inc","textures.inc"] )
	scene.render('regions.png',width=800,height=600 ,antialiasing=.01,remove_temp=False)