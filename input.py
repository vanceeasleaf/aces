#encoding:utf8
import sys
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
	if method=="nvt":
		if(lx/deta/2<upP):exit("upP is too large!")
	if method=="muller":
		if(lx/deta/4<upP):exit("upP is too large!")
	fixl1=xlo-deta;fixl2=fixl1+deta*wfix;
	fixr2=xhi+deta;fixr1=fixr2-deta*wfix;
	fixd1=ylo;fixd2=fixd1+deta;
	fixu2=yhi;fixu1=fixu2-deta;
	hotl1=[0.0]*nstat;hotr2=[0.0]*nstat;
	hotl2=[0.0]*nstat;hotr1=[0.0]*nstat;
	cold=['cold']*nstat;hot=['hot']*nstat
	hotl1[0]=fixl2;hotl2[0]=hotl1[0]+hdeta;
	hotr2[0]=fixr1;hotr1[0]=hotr2[0]-hdeta;
	cold[0]="cold0";
	hot[0]="hot0";

	for i in range(1,nstat):
	    hotl1[i]=hotl2[i-1];hotl2[i]=hotl1[i]+hdeta;
	    hotr2[i]=hotr1[i-1];hotr1[i]=hotr2[i]-hdeta;
	    cold[i]="cold%d"%i;
	    hot[i]="hot%"%i;
	downP=upP;
	down11=xlo+deta*downP;
	down12=down11+deta;
	down22=xhi-deta*downP;
	down21=down22-deta;
	up12=xlo+lx/2-deta*upP;
	up11=up12-deta;
	up21=xlo+lx/2+deta*upP;
	up22=up21+deta;
	lp=up11-down11;



	icoldl=xlo;
	icoldr=icoldl+nswap*deta;
	ihotl=xlo+lx/2;
	ihotr=ihotl+nswap*deta;
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

	#regions and groups
	if method=='nvt':
		print "region	stayl    block  %f %f INF INF INF  INF units box"%(fixl1,fixl2)
		print "region	stayr    block  %f	%f INF INF INF  INF units box"%(fixr1,fixr2)
		for i in range(nstat):
			print "region	%s	block   %f  %f INF INF INF  INF units box"%(cold[i],hotl1[i],hotl2[i])
			print "group	%s	region %s"%(cold[i],cold[i])
			print "region	%s	block	%f	%f	INF INF INF  INF units box"%(hot[i],hotr1[i],hotr2[i])
			print "group	%s	region	%s"%(hot[i],hot[i])


		if fixud:
			print "region	re_nve    block  %f	%f %f	%f INF  INF units box"%(hotl2[-1],hotr1[-1],fixd2,fixu1)
			print "region	stayu	block INF INF %f	%f INF  INF units box"%(fixu1,fixu2)
			print "region	stayd	block INF INF %f	%f INF  INF units box"%(fixd1,fixd2)
			print "region     stay    union  4 stayl stayr stayu stayd"
		else:
			print "region	re_nve    block  %f	%f INF INF INF  INF units box"%(hotl2[-1],hotr1[-1])
			print "region     stay    union  2 stayl stayr "
		print "group            stay    region stay"
		print "group            g_nve    region re_nve"
		s="region	main union %d re_nve "%(nstat*2+1)
		for i in range(nstat):
			s+="%s %s"%(cold[i],hot[i])
		print s
		print "group main region main"
	if method=='muller' or method =='inject':
		print "region            up1    block   %f %f  INF INF INF INF  units box"%(up11,up12)
		print "region            up2    block   %f %f  INF INF INF INF  units box"%(up21,up22)
		print "region            up  union 2 up1 up2"
		print "region            down1  block  %f %f   INF INF INF INF units box"%(down11,down12)
		print "region            down2  block  %f %f   INF INF INF INF units box"%(down21,down22)
		print "region            down union 2 down1 down2"
	print "region	hot block	%f %f   INF INF INF INF units box"%(ihotl,ihotr)
	print "region	cold block %f	%f	INF INF INF INF  units box"%(icoldl,icoldr)
	#computes
	print "compute           ke  all  ke/atom"
	print "compute           pe  all  pe/atom"
	print "compute         stress all stress/atom virial"
	print "compute jflux all heat/flux ke pe stress"
	if method=='muller' or method=='inject':
		print "compute           up_temp    all  temp/region up"
		print "compute           down_temp  all  temp/region down"
	#variables
	if method!='greenkubo':
		print "variable   jx atom (c_pe+c_ke)*vx-(c_stress[1]*vx+c_stress[4]*vy+c_stress[5]*vz)/%f"%nktv
		print "variable   jy atom (c_pe+c_ke)*vy-(c_stress[4]*vx+c_stress[2]*vy+c_stress[6]*vz)/%f"%nktv
		print "variable   jz atom (c_pe+c_ke)*vz-(c_stress[5]*vx+c_stress[6]*vy+c_stress[3]*vz)/%f"%nktv
		print "variable temp atom c_ke/(1.5*%f)"%kb
		print "variable te atom c_ke+c_pe"
		print "variable jcx atom v_te*vx"
		print "variable jcy atom v_te*vy"
		print "variable jcz atom v_te*vz"
		print "variable          delta_temp   equal  c_up_temp-c_down_temp"

	#init atoms to T
	print masses
	print potential
	print "reset_timestep 0"
	print "velocity all create %f %d mom yes rot yes dist gaussian"%(T,seed)
	if method=='nvt':
		print "velocity stay set 0 0 0"
		print "fix getEqu main nvt temp %f %f %f"%(T,T,dtime)
	else:
		print "fix getEqu  all  nvt temp %f %f %f"%(T,T,dtime)
	print "run %d"%equTime
	print "unfix getEqu"
	print "reset_timestep 0"
	if method=='nvt':
		if langevin==0:
			print "fix     nve  g_nve  nve"
		else:
			print "fix     nve  main  nve"
	else:
		if nvt==1:#a key to elimilate the heat of numeric
			print "fix getEqu  all  nvt temp %f %f %f"%(T,T,dtime)
		else:
			print "fix   nve  all  nve"
	print "fix    flux_out  all  ave/time  1  %d  %d  c_jflux[1]  c_jflux[2] c_jflux[3] file  flux.txt "%(aveRate,aveRate)
	if method=='nvt':
		if langevin==1:
			tmstat='langevin';
			r='%d tally yes'%seed;
		else:
			tmstat='nvt temp'
			r=""
		for i in range(nstat):
			print "fix   %s %s %s  %f %f  %s %s"%(hot[i],hot[i],tmstat,Thi,Thi,dtime,r)
			print "fix   %s %s %s  %f %f  %s %s"%(cold[i],cold[i],tmstat,Tlo,Tlo,dtime,r)
			if langevin==0:
				print "compute my%s %s temp/com"%(hot[i],hot[i])
				print "fix_modify %s temp my%s"%(hot[i],hot[i])
				print "compute my%s %s temp"%(cold[i],cold[i])
				print "compute_modify my%s extra 0"%(cold[i])
				print "fix_modify %s temp my%s"%(cold[i],cold[i])
			else:
				print "compute my%s %s temp"%(hot[i],hot[i])
				print "fix_modify %s temp my%s"%(hot[i],hot[i])
				print "compute my%s %s temp"%(cold[i],cold[i])
				print "fix_modify %s temp my%s"%(cold[i],cold[i])
		shot="variable hot equal 0"
		scold="variable cold equal 0"
		for i in range(nstat):
			shot+="+f_%s"%hot[i]
			scold+="+f_%s"%cold[i]
		print shot
		print scold
		print "fix j_hot all ave/time 1 %d  %d v_hot  v_cold file nvtWork.txt "%(aveRate,aveRate)

	if method != "greenkubo":
		print "fix	temp_profile    all    ave/spatial  1  %d  %d  x  lower  %f      v_temp  v_jx file  tempProfile.txt  norm sample units box"%(aveRate,aveRate,deta)
		# 输出热流空间分布,不计算热导率
		if jprofile==1:
			print "dump jprofile all custom %d jprofile.txt id v_jx v_jy v_jz v_temp v_jcx v_jcy v_jcz vx vy vz x y z"%(dumpRate)
			print "dump_modify  jprofile sort id"
	else:
		v=lx*ly*lz
		factor=corRate*timestep/(v*kb*T*T)*zfactor*tcfactor
		#/* 用傅里叶变换热流来计算热流关联函数*/
		if fourierTc==1:
			print "fix j_out all ave/time 1 1 1 c_jflux[1] c_jflux[2] c_jflux[3]  file  jin.txt "
		#/* 用分子模拟论坛的compute扩展来计算热流关联函数,比lammps自带的更精确*/
		if computeTc==1:
			rfactor=tcfactor*zfactor
			if gstart!=20000:gstart=20000
			print "variable          factor_ac equal 1.0"
			print "variable          factor_tc equal %f"%(rfactor)
			print "compute           tc all tc c_thermo_temp c_jflux v_factor_ac v_factor_tc x first %d 0 500000"%(gstart)
			print "fix               tc_out  all  ave/time  1  1  1  c_tc   file  kappa.txt"
		if fourierTc!=1 and computeTc!=1:
			print "fix ss all ave/correlate %d  %d  %d c_jflux[1] c_jflux[2] c_jflux[3] type auto ave running"%(corRate,corNum,aveRate)
			print "variable k11 equal trap(f_ss[3])*%f"%(factor)
			print "variable k22 equal trap(f_ss[4])*%f"%(factor)
			print "variable k33 equal trap(f_ss[5])*%f"%(factor)
			print "fix output all ave/time 1  1 %d v_k11  v_k22  v_k33 file kappa.txt"%(aveRate)
			#/* 定时输出热流自关联函数*/
			if jcf==1:
				print "fix out all ave/time %d  1   %d f_ss[3] f_ss[4] f_ss[5] mode vector file jcf*.txt"%(aveRate,aveRate)
	if method=='muller':
		Nswapbins=2*int(lx/(2*nswap*deta))
		factor=lx/(excRate*timestep)/(2*S)/(1/lp)*zfactor*tcfactor
		print "fix delta_out  all  ave/time  1  %d  %d  v_delta_temp   file  delta_temp.txt"%(aveRate,aveRate)
		print "fix heat_swap   all  thermal/conductivity  %d  x   %d"%(excRate,Nswapbins)
		print "fix e_exchange  all  ave/time  %d  %d  %d  f_heat_swap  file  fileSwap"%(excRate,excNum,aveRate)
		print "variable thermal_conductivity equal f_e_exchange/(1e-10+f_delta_out)*%f"%(factor)
		print "fix thermal_conductivity_out  all  ave/time  %d  1   %d  v_thermal_conductivity   file  thermal_conductivity.txt"%(aveRate,aveRate)
	if method=='inject':
		factor=1/(2*S)/(1/lp)*zfactor*tcfactor
		print "fix delta_out  all  ave/time  1  %d  %d  v_delta_temp   file  delta_temp.txt"%(aveRate,aveRate)
		print "fix               hot   all  heat  %d   %f   region hot"%(excRate,swapEnergyRate)
		print "fix               cold  all  heat  %d  -%f   region cold"%(excRate,swapEnergyRate)
		print "variable          thermal_conductivity equal %f/(1e-10+f_delta_out)*%f"%(swapEnergyRate,factor)
		print "fix thermal_conductivity_out  all  ave/time  %d  1   %d  v_thermal_conductivity   file  thermal_conductivity.txt"%(aveRate,aveRate)


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
