#encoding:utf8
import sys,json
class nvtDevice:
	def __init__(self,hook,m):
		self.para=[m.box,m.deta,m.wfix,m.nstat,m.upP,m.hdeta,m.fixud,m.langevin]
		xlo,xhi,ylo,yhi,zlo,zhi,lx,ly,lz=m.box
		if(lx/m.deta/2<m.upP):
			print "upP is too large!"
			sys.exit()
		m.Thi=m.T+m.dT
		m.Tlo=m.T-m.dT
		self.__dict__=dict(self.__dict__,**m.__dict__)
		self.m=m
		self.regions=[]
		hook.addAction('region',self.renderRegion)
		hook.addAction('variable',self.renderVariable)
		hook.addAction('equ',self.renderEqu)
		hook.addAction('elimination',self.renderElim)
		hook.addAction('temp',self.renderTemp)
		hook.addAction('flux',self.renderFlux)
		 
	def addBox(self,name,paras):
		p1,p2=paras		
		print "region	%s	block   %s  %s %s  %s %s  %s units box"%(name,p1[0],p2[0],p1[1],p2[1],p1[2],p2[2])
		self.regions.append({'name':name,'type':'box','dim':paras})
	def renderRegion(self):
		box,deta,wfix,nstat,upP,hdeta,fixud,langevin=self.para
		xlo,xhi,ylo,yhi,zlo,zhi,lx,ly,lz=box
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
		INF='INF'
		self.addBox('stayl',([fixl1,INF,INF],[fixl2,INF,INF]))
		self.addBox('stayr',([fixr1,INF,INF],[fixr2,INF,INF]))
		for i in range(nstat):
			self.addBox(cold[i],([hotl1[i],INF,INF],[hotl2[i],INF,INF]))
			print "group	%s	region %s"%(cold[i],cold[i])
			self.addBox(hot[i],([hotr1[i],INF,INF],[hotr2[i],INF,INF]))
			print "group	%s	region	%s"%(hot[i],hot[i])

		self.hot=hot
		self.cold=cold
		self.hotl2=hotl2
		self.hotr1=hotr1
		if fixud:
			self.addBox('re_nve',([hotl2[-1],fixd2,INF],[hotr1[-1],fixu1,INF]))
			self.addBox('stayu',([INF,fixu1,INF],[INF,fixu2,INF]))
			self.addBox('stayd',([INF,fixd1,INF],[INF,fixd2,INF]))
			print "region     stay    union  4 stayl stayr stayu stayd"
		else:
			self.addBox('re_nve',([hotl2[-1],INF,INF],[hotr1[-1],INF,INF]))
			print "region     stay    union  2 stayl stayr "
		print "group            stay    region stay"
		print "group            g_nve    region re_nve"
		s="region	main union %d re_nve "%(nstat*2+1)
		for i in range(nstat):
			s+="%s %s"%(cold[i],hot[i])
		print s
		print "group main region main"
		f=open('regions.txt','w')
		f.write(json.dumps(self.regions))
	def renderVariable(self):
		print "variable   jx atom (c_pe+c_ke)*vx-(c_stress[1]*vx+c_stress[4]*vy+c_stress[5]*vz)/%f"%self.nktv
		print "variable   jy atom (c_pe+c_ke)*vy-(c_stress[4]*vx+c_stress[2]*vy+c_stress[6]*vz)/%f"%self.nktv
		print "variable   jz atom (c_pe+c_ke)*vz-(c_stress[5]*vx+c_stress[6]*vy+c_stress[3]*vz)/%f"%self.nktv
		print "variable temp atom c_ke/(1.5*%f)"%self.kb
		print "variable te atom c_ke+c_pe"
		print "variable jcx atom v_te*vx"
		print "variable jcy atom v_te*vy"
		print "variable jcz atom v_te*vz"
		print "variable          delta_temp   equal  c_up_temp-c_down_temp"
	def renderEqu(self):
		#print "velocityRamp g_nve create %f %d  vx %f %f x %f %f mom yes rot yes dist gaussian"%(self.T,self.seed,self.Tlo,self.Thi,self.hotl2[-1],self.hotr1[-1])
		#print "velocity g_nve create %f %d   mom yes rot yes dist gaussian"%(self.T,self.seed)
		print "velocity stay set 0 0 0"
		print "fix getEqu main nvt temp %f %f %f"%(self.T,self.T,self.dtime)

	def renderElim(self):
		bath="nve"
		if self.m.nvt:bath="nvt temp %f %f %f"%(self.T,self.T,self.dtime)
		if self.langevin==0:
			print "fix     nve  g_nve  %s"%bath
		else:
			print "fix     nve  main  %s"%bath
	def renderTemp(self):
		hot=self.hot
		cold=self.cold
		if self.langevin==1:
			tmstat='langevin';
			r='%d tally yes'%self.seed;
		else:
			tmstat='nvt temp'
			r=""
		for i in range(self.nstat):
			
			if self.langevin==0:
				print "fix   %s %s %s  %f %f  %s %s"%(hot[i],hot[i],tmstat,self.Thi,self.Thi,self.dtime,r)
				print "fix   %s %s %s  %f %f  %s %s"%(cold[i],cold[i],tmstat,self.Tlo,self.Tlo,self.dtime,r)
				print "compute my%s %s temp/com"%(hot[i],hot[i])
				print "fix_modify %s temp my%s"%(hot[i],hot[i])
				print "compute my%s %s temp"%(cold[i],cold[i])
				print "compute_modify my%s extra 0"%(cold[i])
				print "fix_modify %s temp my%s"%(cold[i],cold[i])
			else:
				print "fix   %s %s %s  %f %f  %s %s"%(hot[i],hot[i],tmstat,self.Thi*3.0/self.dimension,self.Thi*3.0/self.dimension,self.dtime,r)
				print "fix   %s %s %s  %f %f  %s %s"%(cold[i],cold[i],tmstat,self.Tlo*3.0/self.dimension,self.Tlo*3.0/self.dimension,self.dtime,r)
				print "compute my%s %s temp"%(hot[i],hot[i])
				print "fix_modify %s temp my%s"%(hot[i],hot[i])
				print "compute my%s %s temp"%(cold[i],cold[i])
				print "fix_modify %s temp my%s"%(cold[i],cold[i])
		shot="variable hot equal 0"
		scold="variable cold equal 0"
		for i in range(self.nstat):
			shot+="+f_%s"%hot[i]
			scold+="+f_%s"%cold[i]
		print shot
		print scold
		print "fix j_hot all ave/time 1 %d  %d v_hot  v_cold file nvtWork.txt "%(self.aveRate,self.aveRate)
	def renderFlux(self):
		print "fix	temp_profile    main    ave/spatial  1  %d  %d  x  lower  %f      v_temp  v_jx file  tempProfile.txt  norm sample units box"%(self.aveRate,self.aveRate,self.deta)
		# 输出热流空间分布,不计算热导率
		if self.jprofile==1:
			print "dump jprofile main custom %d jprofile.txt id v_jx v_jy v_jz v_temp v_jcx v_jcy v_jcz vx vy vz x y z"%(self.dumpRate)
			print "dump_modify  jprofile sort id"

