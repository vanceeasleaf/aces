#encoding:utf8
from math import sin,cos,atan,pi,sqrt
from os.path import *
import os,imp
root=abspath(dirname(realpath(__file__))+'/../../');
default= imp.load_source('default', root+'/default.py') 
unitcell= imp.load_source('unitcell', root+'/UnitCell/unitcell.py') 
from ase.io.vasp import write_vasp
UnitCell=unitcell.UnitCell
class structure:
	def __init__(self,home,opt):
		self.home=home
		self.__dict__=dict(self.__dict__,**default.default)# all the values needed
		self.potential='pair_style        tersoff\npair_coeff      * * %s/potentials/BNC.tersoff  C N'%home
		self.masses=""
		self.dump="dump_modify dump1 element C N"
		self.latx=11;self.laty=1;self.latz=1;
		self.__dict__=dict(self.__dict__,**opt)
	def structure(self):
		home=self.home
		latx,laty,latz,bond=[self.latx,self.laty,self.latz,self.bond]
		pos1=[
		6.928400,13.000369,0.000000
		,7.794450,16.500469,0.000000
		,9.526550,10.500299,0.000000
		,11.258650,17.500498,0.000000
		]
		phi=pi/2-atan((pos1[4]-pos1[1])/(pos1[3]-pos1[0]));
		bond=sqrt((pos1[4]-pos1[1])*(pos1[4]-pos1[1])+(pos1[3]-pos1[0])*(pos1[3]-pos1[0]))*1.42;
		dx=sqrt(3)*bond;
		dy=3*bond;
		pos2=[
		7.794450,13.500383,0.000000
		,6.928400,12.000340,0.000000
		,12.124700,13.000369,0.000000
		,11.258650,13.500383,0.000000
		,10.392600,13.000369,0.000000
		,10.392600,12.000340,0.000000
		,8.660500,13.000369,0.000000
		,11.258650,14.500412,0.000000
		,12.124700,12.000340,0.000000
		,9.526550,14.500412,0.000000
		,9.526550,13.500383,0.000000
		,8.660500,12.000340,0.000000
		,7.794450,14.500412,0.000000
		,12.990750,13.500383,0.000000
		,12.990750,14.500412,0.000000
		,9.526550,11.500326,0.000000
		,7.794450,11.500326,0.000000
		,12.124700,10.000284,0.000000
		,11.258650,11.500326,0.000000
		,11.258650,10.500299,0.000000
		,10.392600,21.000597,0.000000
		,14.722850,19.500553,0.000000
		,11.258650,20.500582,0.000000
		,12.124700,18.000511,0.000000
		,11.258650,19.500553,0.000000
		,9.526550,19.500553,0.000000
		,9.526550,20.500582,0.000000
		,8.660500,18.000511,0.000000
		,8.660500,19.000540,0.000000
		,12.990750,19.500553,0.000000
		,10.392600,18.000511,0.000000
		,10.392600,19.000540,0.000000
		,12.990750,20.500582,0.000000
		,13.856800,18.000511,0.000000
		,12.124700,19.000540,0.000000
		,13.856800,19.000540,0.000000
		,13.856800,16.000454,0.000000
		,12.990750,17.500498,0.000000
		,12.124700,15.000426,0.000000
		,10.392600,16.000454,0.000000
		,11.258650,16.500469,0.000000
		,12.124700,16.000454,0.000000
		,12.990750,16.500469,0.000000
		,8.660500,15.000426,0.000000
		,8.660500,16.000454,0.000000
		,10.392600,15.000426,0.000000
		,9.526550,16.500469,0.000000
		,9.526550,17.500498,0.000000
		,6.928400,13.000369,0.000000
		,7.794450,16.500469,0.000000
		,9.526550,10.500299,0.000000
		,11.258650,17.500498,0.000000
		]
		rpos=[0.0]*52
		for i in range(52):
			rpos[i]=self.rotate(phi,pos2[i*3],pos2[i*3+1]);
		minx=100000;
		miny=100000;
		for i in range(52):
			minx=min(rpos[i][0],minx);
			miny=min(rpos[i][1],miny);
		for i in range(52):
			rpos[i][0]-=minx;
			rpos[i][1]-=miny;
		file=open("cell.in","w");
		content="""1
%d 0 0 
0 %d 0
0 0 100
51 1
"""%(dx,dy)
		for i in range(52):
			content+="%f\t%f\t%f\n"%(rpos[i][0]*1.42/dx,rpos[i][1]*1.42/dy,0.000000)
		
		file.write(content)
		content="""7
y
cell.in
y
%d %d %d
1
3
C N
0
CN.xyz
structure
map.in
"""%(latx,laty,latz)
		file.close()		
		file=open("in.disp","w");
		file.write(content)
		file.close()
		os.popen(home+"/latgen <in.disp").read();
		print 'read_data structure'
		
	def rotate(self,phi,x,y):
		x1=cos(phi)*x-sin(phi)*y
		y1=sin(phi)*x+cos(phi)*y
		return [x1,y1]
	def write(self):
		self.atoms.write("CN.xyz")
		write_vasp("POSCAR",self.atoms,sort="True",direct=True,vasp5=True)
		poscar = open("POSCAR")
		unit_cell = UnitCell(poscar)
		unit_cell.num_atom_types=2
		lammps=open("structure","w")
		lammps.write(unit_cell.output_lammps())
		lammps.close()