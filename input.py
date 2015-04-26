#encoding:utf8
import sys
from devices import nvtDevice,mpDevice,ijDevice,gkDevice
from aces.Units import Units
from aces.tools import *
import aces.config as config
from ase.io import read
from ase.io.vasp import write_vasp
def getboxrange():
	file=open("minimize/range");
	for i in range(5):
		file.next()
	xlo,xhi=map(float,file.next().split()[:2])
	ylo,yhi=map(float,file.next().split()[:2])
	zlo,zhi=map(float,file.next().split()[:2])
	return (xlo,xhi,ylo,yhi,zlo,zhi);

def getxrange():
    file=open('minimize/minimize.xyz');
    n=file.next().split()[0];n=int(n)
    file.next()
    xmin=100000;xmax=-100000;
    ymin=100000;ymax=-100000;
    zmin=100000;zmax=-100000;
    for i in range(n):
        label,x,y,z=file.next().split();
        x,y,z=map(float,[x,y,z]);
        xmin=min(x,xmin);xmax=max(x,xmax);
        ymin=min(y,ymin);ymax=max(y,ymax);
        zmin=min(z,zmin);zmax=max(z,zmax); 
    return (xmin,xmax,ymin,ymax,zmin,zmax);
    
def postMini(xp,yp,zp,enforceThick,thick):
	xlo,xhi,ylo,yhi,zlo,zhi=getboxrange();
	xlo0,xhi0,ylo0,yhi0,zlo0,zhi0=getxrange();
	if(xp==0):
		xlo=xlo0;xhi=xhi0;
	if(yp==0):
		ylo=ylo0;yhi=yhi0;
	if(zp==0):
		zlo=zlo0;zhi=zhi0;
	lx=xhi-xlo;ly=yhi-ylo;lz=zhi-zlo;
	if(enforceThick):zfactor=lz/thick;
	else:zfactor=1;
	S=ly*lz;
	return lx,ly,lz,zfactor,S,xlo,xhi,ylo,yhi,zlo,zhi



class Hook:
	def __init__(self):
		self.labels={}
	def addAction(self,label,function):
		if not self.labels.has_key(label):
			self.labels[label]=[]
		self.labels[label].append(function)
	def doAction(self,label):
		if not self.labels.has_key(label):return
		for key in self.labels[label]:
			key()
def mdTc(m):
	units=Units(m.units)
	m.kb=units.boltz
	m.nktv=units.nktv2p
	if(m.method=="nvt"):m.xp=0;
	lx,ly,lz,m.zfactor,m.S,xlo,xhi,ylo,yhi,zlo,zhi=postMini(m.xp,m.yp,m.zp,m.enforceThick,m.thick)
	m.dtime=m.timestep*100;
	m.tcfactor=units.tcfactor;
	m.excNum=m.aveRate/m.excRate;
	m.swapEnergyRate=m.swapEnergy/(m.excRate*m.timestep);

	m.box=(xlo,xhi,ylo,yhi,zlo,zhi,lx,ly,lz)
	hook=Hook()
	if m.method=='nvt':
		device=nvtDevice(hook,m)
	elif m.method=='muller':
		device=mpDevice(hook,m)
	elif m.method=='inject':
		device=ijDevice(hook,m)
	elif m.method=='greenkubo':
		m.corNum=m.aveRate/m.corRate;
		device=gkDevice(hook,m)

	#settings
	print "units %s"%units
	print "dimension 3"
	pbcx=pbcy=pbcz='s'
	if m.xp==1:pbcx='p'
	if m.yp==1:pbcy='p'
	if m.zp==1:pbcz='p'
	print "boundary %s %s %s"%(pbcx,pbcy,pbcz)
	print "atom_style atomic"
	print "read_restart   minimize/restart.minimize"
	print "change_box	all	boundary %s %s %s"%(pbcx,pbcy,pbcz)
	print "lattice fcc 5" #needed to define the regions
	print "thermo %d"%m.dumpRate
	print "thermo_modify     lost warn"
	print "timestep %f"%m.timestep
	#regions and groups
	hook.doAction('region')
	#computes
	print "compute           ke  all  ke/atom"
	print "compute           pe  all  pe/atom"
	print "compute         stress all stress/atom virial"
	print "compute jflux all heat/flux ke pe stress"
	hook.doAction('compute')
	#variables
	hook.doAction('variable')
	#init atoms to T
	print m.masses
	print m.potential
	print "reset_timestep 0"
	print "velocity all create %f %d mom yes rot yes dist gaussian"%(m.T,m.seed)
	hook.doAction('equ')
	print "run %d"%m.equTime
	print "unfix getEqu"
	print "reset_timestep 0"
	hook.doAction('elimination')
	print "fix    flux_out  all  ave/time  1  %d  %d  c_jflux[1]  c_jflux[2] c_jflux[3] file  flux.txt "%(m.aveRate,m.aveRate)
	hook.doAction('temp')
	hook.doAction('flux')
	hook.doAction('swap')


	#/* 定时输出dump文件并按id排序*/
	if(m.dumpxyz):
		print "dump dump1 all atom %d dump.lammpstrj"%(m.dumpRate)
		print "dump_modify  dump1 sort id"


	#/* 定时输出速度文件用于计算速度关联函数*/
	if(m.dumpv):
		print "dump dump2 all custom %d dump.velocity type vx vy vz"%(m.dumpRate)
		print "dump_modify  dump2 sort id"

	print "run	%d"%(m.runTime)
def bte(m):
	coordination=phontsAtoms()
	content0="species %d\n"%(len(m.elements))+m.phontsmasses+"""
D3_cutoff 9.0
delta 0.0005
numerical_2der T

iter_steps 10

AbInitio  T F 
FP_interface LAMMPS
phonons_only T
Lattice  1.0
%s
end
"""%coordination
	write(content0,'phonons_input.dat')
	passthru(config.phonts) # generate many displacement files
	mkdir('lammps');cd('lammps')
	content="units %s\n"%m.units
	content+="""atom_style      charge
dimension       3
boundary        p p p 
read_data       GENERIC
%s
%s 
neighbor        1.1 bin
neigh_modify    every 1 delay 1 check yes
dump 1 all custom 1 *.dump id  fx fy fz 
dump_modify 1 format "%%d %%30.20f %%30.20f %%30.20f"
dump_modify  1 sort id
run 0
"""%(m.masses,m.potential)
	shell_exec("mv ../*.str .")
	strs=shell_exec("ls *.str").split("\n")
	for str in strs:
		dir=str.replace("str","dir")
		mkdir(dir)
		write(content.replace("GENERIC",str),dir+"/in")
		mv(str,"%s/%s"%(dir,str))
		cd(dir)
		passthru(config.lammps+" <in >out.dat")
		cd('..')
	dirs=shell_exec("ls |grep dir").split("\n")
	for dir in dirs:
		print dir
		cp(dir+"/0.dump","../"+dir.replace("dir","out"))
	cp("1.0000.dir/out.dat","../1.0000.out")
	cd('..')
	content0=content0.replace('AbInitio  T F','AbInitio  F T')
	write(content0,'phonons_input.dat')
	passthru(config.phonts)
def input(m):
	if m.bte==True:bte(m)
	elif m.bte=="phonopy":phopy(m)
	elif m.bte=="shengbte":shengbte(m)
	elif m.bte=="correlation":correlation(m)
	else:mdTc(m)
def phopy(m):
	# from minimized structure generate POSCAR
	atoms=read('minimize/range',format='lammps')
	s=atoms.numbers
	symbols=[m.elements[i-1] for i in s ]
	atoms.set_chemical_symbols(symbols)
	write_vasp("POSCAR",atoms,sort="True",direct=True,vasp5=True)
	#generate supercells
	dim=' '.join(str(i) for i in m.supercell)
	passthru(config.phonopy+"-d --dim='%s'"%(dim))
	files=shell_exec("ls *-*").split('\n')
	content="units %s\n"%m.units
	content+="""atom_style      atomic
dimension       3
boundary        p p p 
read_data       structure
%s
%s 
neighbor        2 bin
neigh_modify    every 1 delay 1 check yes
dump 1 all custom 1 dump.force id  fx fy fz xs ys zs
dump_modify 1 format "%%d %%30.20f %%30.20f %%30.20f %%30.20f %%30.20f %%30.20f"
dump_modify  1 sort id
run 0
"""%(m.masses,m.potential)
	from aces.UnitCell.unitcell import UnitCell
	cmd=config.phonopy+"-f "
	maindir=shell_exec('pwd')
	for file in files:
		print file
		dir="dirs/dir_"+file
		mkdir(dir)
		#POSCAR
		mv(file,dir)
		cd(dir)
		#generate structure
		poscar = open(file)
		unit_cell = UnitCell(poscar)
		unit_cell.num_atom_types=len(m.elements)
		write(unit_cell.output_lammps(),"structure")
		#generate in
		write(content,"in")
		#generate dump.force
		shell_exec(config.lammps+" < in >log.out")
		#generate vasprun.xml
		f=open('dump.force')
		for i in range(9):f.next()
		forces=""
		poses=""
		for line in f:
			line=line.split()
			forces+="<v>  %s %s %s </v>\n"%tuple(line[1:4])
			poses+="<v>  %s %s %s  </v>\n"%tuple(line[4:8])
		vasprun='<root><calculation><varray name="forces" >\n'
		vasprun+=forces
		vasprun+='</varray>\n<structure><varray name="positions">\n'+poses
		vasprun+='</varray></structure></calculation></root>\n'
		write(vasprun,'vasprun.xml')
		cmd+=dir+'/vasprun.xml '
		cd(maindir)
	#generate FORCE_SETS
	passthru(cmd)
	#generate mesh.conf
	mesh="""DIM = %s
ATOM_NAME = %s
MP = %s
EIGENVECTORS=.TRUE.
FORCE_CONSTANTS = WRITE
"""%(dim,' '.join(m.elements),' '.join(str(i) for i in m.kpoints))
	write(mesh,'mesh.conf')
	import matplotlib
	matplotlib.use('Agg')
	from matplotlib import pyplot as pl
	passthru(config.phonopy+" --dos  mesh.conf")

	import numpy as np
	xx=np.loadtxt('partial_dos.dat',skiprows=1)
	ndos=len(line)-1
	freq=xx[:,0]
	pdos=xx[:,1:]
	ndos=len(pdos[0,:])
	pl.figure()
	for i in range(ndos):
		pl.plot(freq,pdos[:,i])
	pl.xlabel('Frequency (THz)')
	pl.ylabel('Partial Density of States')
	pl.savefig('partial_dos.png',bbox_inches='tight',transparent=True) 
	pl.figure()
	pl.plot(freq,np.sum(pdos,axis=1),color='red')
	pl.xlabel('Frequency (THz)')
	pl.ylabel('Density of States')
	pl.savefig('total_dos.png',bbox_inches='tight',transparent=True) 
	#calculate paticipation ratio
	from aces.binary import pr
	pr()
	#plot
	xs=[];ys=[]
	for line in open('pr.txt'):
		x,y=map(float,line.split())
		xs.append(x);ys.append(y)
	write("%s"%(sum(ys)/len(ys)),"ave_pr.txt")
	pl.figure()
	pl.plot(xs,ys,'.',color='red')
	pl.ylim([0.0,1.0])
	pl.xlabel('Frequency (THz)')
	pl.ylabel('Paticipation Ratio')
	pl.savefig('Paticipation_atio.png',bbox_inches='tight',transparent=True) 
def shengbte(m):
	pass	
def correlation(m):
	units=Units(m.units)
	m.kb=units.boltz
	m.nktv=units.nktv2p
	if(m.method=="nvt"):m.xp=0;
	lx,ly,lz,m.zfactor,m.S,xlo,xhi,ylo,yhi,zlo,zhi=postMini(m.xp,m.yp,m.zp,m.enforceThick,m.thick)
	m.dtime=m.timestep*100;
	m.tcfactor=units.tcfactor;
	m.excNum=m.aveRate/m.excRate;
	m.swapEnergyRate=m.swapEnergy/(m.excRate*m.timestep);
	f=open("correlation.lmp","w")
	print >>f,"units %s"%units
	print >>f,"dimension 3"
	pbcx=pbcy=pbcz='s'
	if m.xp==1:pbcx='p'
	if m.yp==1:pbcy='p'
	if m.zp==1:pbcz='p'
	print >>f,"boundary %s %s %s"%(pbcx,pbcy,pbcz)
	print >>f,"atom_style atomic"
	print >>f,"read_restart   minimize/restart.minimize"
	print >>f,"change_box	all	boundary %s %s %s"%(pbcx,pbcy,pbcz)
	print >>f,"lattice fcc 5" #needed to define the regions
	print >>f,"thermo %d"%m.dumpRate
	print >>f,"thermo_modify     lost warn"
	print >>f,m.masses
	print >>f,m.potential
	print >>f,"timestep %f"%m.timestep
	print >>f,"reset_timestep 0"
	print >>f,"velocity all create %f %d mom yes rot yes dist gaussian"%(m.T,m.seed)
	print >>f,"fix getEqu  all  nvt temp %f %f %f"%(m.T,m.T,m.dtime)
	print >>f,"run %d"%m.equTime
	print >>f,"unfix getEqu"
	print >>f,"reset_timestep 0"
	print >>f,"fix nve all nve"
	print >>f,"dump lala all custom %s velocity.txt id type vx vy vz"%m.Cinterval
	print >>f,"dump_modify  lala sort id"
	print >>f,"run %s"%m.Ctime
	f.close()
	passthru(config.lammps+" <correlation.lmp >out.dat")
	from aces.vdos import vdos
	vdos(m.timestep)
	rm("velocity.txt")
def phontsAtoms():
	atoms=read('minimize/range',format='lammps')
	cell=atoms.get_cell()
	content="cell %f %f %f\n"%(cell[0][0],cell[1][1],cell[2][2])
	content+="natoms %d\n"%(len(atoms))
	content+="fractional\n"
	pos=atoms.get_scaled_positions()
	for i,atom in enumerate(atoms):
		if atom.symbol=='H':label='C'
		if atom.symbol=='He':label='N'
		content+="%s %s\n"%(label,' '.join(["%s"%x for x in pos[i]]))
	return content		
if __name__=='__main__':
	units ,xp ,yp ,zp ,dumpRate ,timestep ,method ,kb ,nktv ,masses,potential ,T ,seed ,dtime ,equTime ,langevin ,nvt ,aveRate ,deta ,jprofile ,dumpRate ,corRate ,computeTc  ,fourierTc ,tcfactor ,gstart ,jcf  ,nswap ,excRate  ,excNum ,swapEnergyRate ,dumpxyz ,dumpv ,runTime,upP,wfix,nstat,enforceThick,thick,Thi,Tlo,hdeta,fixud=sys.argv[1:]
	#input(units ,xp ,yp ,zp ,dumpRate ,timestep ,method ,kb ,nktv ,masses,potential ,T ,seed ,dtime ,equTime ,langevin ,nvt ,aveRate ,deta ,jprofile ,corRate ,computeTc  ,fourierTc ,tcfactor ,gstart ,jcf  ,nswap ,excRate  ,excNum ,swapEnergyRate ,dumpxyz ,dumpv ,runTime,upP,wfix,nstat,enforceThick,thick,Thi,Tlo,hdeta,fixud)
