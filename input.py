#encoding:utf8
import sys
from nvtDevice import nvtDevice
from mpDevice import mpDevice,ijDevice
from gkDevice import gkDevice
def exit(info):
	print info
	sys.exit();
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





def input(units ,xp ,yp ,zp ,dumpRate ,timestep ,method ,kb ,nktv ,masses,potential ,T ,seed ,dtime ,equTime ,langevin ,nvt ,aveRate ,deta ,jprofile,corRate ,computeTc  ,fourierTc ,tcfactor ,gstart ,jcf  ,nswap ,excRate  ,excNum ,swapEnergyRate ,dumpxyz ,dumpv ,runTime,upP,wfix,nstat,enforceThick,thick,Thi,Tlo,hdeta,fixud):
	xp ,yp ,zp ,dumpRate,seed,equTime,langevin ,nvt ,aveRate ,jprofile,corRate ,computeTc,fourierTc,gstart ,jcf  ,nswap ,excRate,excNum,dumpxyz ,dumpv ,runTime,upP,wfix,nstat,enforceThick,fixud=map(int,[xp ,yp ,zp ,dumpRate,seed,equTime,langevin ,nvt ,aveRate ,jprofile  ,corRate ,computeTc,fourierTc,gstart ,jcf  ,nswap ,excRate,excNum,dumpxyz ,dumpv ,runTime,upP,wfix,nstat,enforceThick,fixud])
	timestep,kb ,nktv,T,dtime,deta,tcfactor ,swapEnergyRate,thick,Thi,Tlo,hdeta=map(float,[timestep,kb ,nktv,T,dtime,deta,tcfactor,swapEnergyRate,thick,Thi,Tlo,hdeta])

	lx,ly,lz,zfactor,S,xlo,xhi,ylo,yhi,zlo,zhi=postMini(xp,yp,zp,enforceThick,thick)
	box=(xlo,xhi,ylo,yhi,zlo,zhi,lx,ly,lz)
	if method=='nvt':
		device=nvtDevice(box,deta,wfix,nstat,upP,hdeta,fixud,langevin,nktv,kb,T,dtime,Thi,Tlo,aveRate,jprofile,dumpRate)
	elif method=='muller':
		device=mpDevice(box,upP,nktv,kb,aveRate,swapEnergyRate,S,zfactor,tcfactor,deta,excRate,excNum,jprofile,dumpRate,timestep,T,dtime,nswap)
	elif method=='inject':
		device=ijDevice(self,box,upP,nktv,kb,aveRate,swapEnergyRate,S,zfactor,tcfactor,deta,excRate,excNum,jprofile,dumpRate,timestep,T,dtime,nswap,swapEnergyRate)
	elif method=='greenkubo':
		corNum=aveRate/corRate;
		device=gkDevice(box,corRate,timestep,kb,T,zfactor,tcfactor,fourierTc,computeTc,gstart,corNum,aveRate,jcf,dtime)

	#settings
	print "units %s"%units
	print "dimension 3"
	pbcx=pbcy=pbcz='s'
	if xp==1:pbcx='p'
	if yp==1:pbcy='p'
	if zp==1:pbcz='p'
	print "boundary %s %s %s"%(pbcx,pbcy,pbcz)
	print "atom_style atomic"
	print "read_restart   minimize/restart.minimize"
	print "change_box	all	boundary %s %s %s"%(pbcx,pbcy,pbcz)
	print "lattice fcc 5" #needed to define the regions
	print "thermo %d"%dumpRate
	print "thermo_modify     lost warn"
	print "timestep %f"%timestep
	device.renderRegion()
	#regions and groups
	#computes
	print "compute           ke  all  ke/atom"
	print "compute           pe  all  pe/atom"
	print "compute         stress all stress/atom virial"
	print "compute jflux all heat/flux ke pe stress"
	device.renderCompute()
	#variables
	device.renderVariable()
	#init atoms to T
	print masses
	print potential
	print "reset_timestep 0"
	print "velocity all create %f %d mom yes rot yes dist gaussian"%(T,seed)
	device.renderEqu()
	print "run %d"%equTime
	print "unfix getEqu"
	print "reset_timestep 0"
	device.renderElim()
	print "fix    flux_out  all  ave/time  1  %d  %d  c_jflux[1]  c_jflux[2] c_jflux[3] file  flux.txt "%(aveRate,aveRate)
	device.renderTemp()
	device.renderFlux()
	device.renderSwap()


	#/* 定时输出dump文件并按id排序*/
	if(dumpxyz):
		print "dump dump1 all atom %d dump.lammpstrj"%(dumpRate)
		print "dump_modify  dump1 sort id"


	#/* 定时输出速度文件用于计算速度关联函数*/
	if(dumpv):
		print "dump dump2 all custom %d dump.velocity type vx vy vz"%(dumpRate)
		print "dump_modify  dump2 sort id"

	print "run	%d"%(runTime)


if __name__=='__main__':
	units ,xp ,yp ,zp ,dumpRate ,timestep ,method ,kb ,nktv ,masses,potential ,T ,seed ,dtime ,equTime ,langevin ,nvt ,aveRate ,deta ,jprofile ,dumpRate ,corRate ,computeTc  ,fourierTc ,tcfactor ,gstart ,jcf  ,nswap ,excRate  ,excNum ,swapEnergyRate ,dumpxyz ,dumpv ,runTime,upP,wfix,nstat,enforceThick,thick,Thi,Tlo,hdeta,fixud=sys.argv[1:]
	input(units ,xp ,yp ,zp ,dumpRate ,timestep ,method ,kb ,nktv ,masses,potential ,T ,seed ,dtime ,equTime ,langevin ,nvt ,aveRate ,deta ,jprofile ,corRate ,computeTc  ,fourierTc ,tcfactor ,gstart ,jcf  ,nswap ,excRate  ,excNum ,swapEnergyRate ,dumpxyz ,dumpv ,runTime,upP,wfix,nstat,enforceThick,thick,Thi,Tlo,hdeta,fixud)
