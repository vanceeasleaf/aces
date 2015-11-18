#this script is used to parrallelize the procedure of getting vasprun.xml
from aces.tools import *
from aces import config
from aces.App import App
import time
import numpy as np
from aces.f import RotateVector

def get_lammps_script(m):
	content="units %s\n"%m.units
	content+="""%s
dimension       3
boundary        p p p 
read_data       structure
%s
%s 
#neighbor        14 bin
#neigh_modify    every 1 delay 1 check yes
dump 1 all custom 1 dump.force id  fx fy fz xs ys zs
dump_modify 1 format "%%d %%f %%f %%f %%f %%f %%f"
dump_modify  1 sort id
run 0
"""%(m.getatomicstyle(),m.masses,m.potential)
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
def run(file):
	
	a=time.time()
	print file
	if not exists(file):
		exit("File not exist!")
	dir="dirs/dir_"+file
	m=App().m #before cd

	if not exists(dir):mkdir(dir)
	mv(file,dir+'/POSCAR')
	#debug(time.strftime('%Y-%m-%d %H:%M:%S'))
	cd(dir)		
	debug('one vasprun mv etc.:%f s'%(time.time()-a))
	a=time.time()
	getVaspRun_lammps(m)
	debug('one vasprun lmp:%f s'%(time.time()-a))
def exe(file):

	passthru(config.python+__file__+" "+file)
if __name__=="__main__":
	run(sys.argv[1])
