#this script is used to parrallelize the procedure of getting vasprun.xml
from aces.tools import *
from aces import config
from aces.App import App
import time
import numpy as np
from aces.f import RotateVector

def get_lammps_script(m):
	content="units %s\n"%m.units
	pbcx=pbcy=pbcz='s'
	if m.xp==1:pbcx='p'
	if m.yp==1:pbcy='p'
	if m.zp==1:pbcz='p'
	b="boundary %s %s %s"%(pbcx,pbcy,pbcz)
	content+="""%s
dimension       3
%s
read_data       structure
%s
%s 
neighbor        2 nsq
neigh_modify    every 1 delay 1 check yes
dump 1 all custom 1 dump.force id  fx fy fz xs ys zs
dump_modify 1 format "%%d %%f %%f %%f %%f %%f %%f"
dump_modify  1 sort id
run 0
"""%(b,m.getatomicstyle(),m.masses,m.potential)
	return content
def getVaspRun_lammps(m):
	#generate structure
	rot=m.POSCAR2data()
	#generate in
	content=get_lammps_script(m)
	write(content,"in")
	#generate dump.force
	shell_exec(config.lammps+" < in >log.out")
	d,p,d1,p1=rot

	#generate vasprun.xml
	f=open('dump.force')
	for i in range(9):f.next()
	forces=""
	poses=""
	for line in f:
		line=line.split()
		force=np.array(map(float,line[1:4]))
		pos=np.array(map(float,line[4:8]))
		force=RotateVector(force,d1,-p1)
		force=RotateVector(force,d,-p)
		
		forces+="<v>  %f %f %f </v>\n"%tuple(force)
		poses+="<v>  %f %f %f  </v>\n"%tuple(pos)
	vasprun='<root><calculation><varray name="forces" >\n'
	vasprun+=forces
	vasprun+='</varray>\n<structure><varray name="positions">\n'+poses
	vasprun+='</varray></structure></calculation></root>\n'
	write(vasprun,'vasprun.xml')
	f.close()
def run():
	a=time.time()
	m=App().m #before cd
	debug('one vasprun mv etc.:%f s'%(time.time()-a))
	a=time.time()
	getVaspRun_lammps(m)
	debug('one vasprun lmp:%f s'%(time.time()-a))
def exe():

	passthru(config.python+__file__)
if __name__=="__main__":
	run()
