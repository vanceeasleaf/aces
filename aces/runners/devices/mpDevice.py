#encoding:utf8
import sys
class mpDevice:
	def __init__(self,hook,m):
		self.para=[m.box,m.upP,m.deta,m.nswap]
		self.__dict__=dict(self.__dict__,**m.__dict__)
		self.lx=m.box[6]
		self.xlo=m.box[0]
		if(self.lx/m.deta/2<m.upP):
			print "upP is too large!"
			sys.exit()
		hook.addAction('region',self.renderRegion)
		hook.addAction('compute',self.renderCompute)
		hook.addAction('variable',self.renderVariable)
		hook.addAction('equ',self.renderEqu)
		hook.addAction('elimination',self.renderElim)
		hook.addAction('flux',self.renderFlux)
		hook.addAction('swap',self.renderSwap)
		self.m=m
	def renderRegion(self):
		box,upP,deta,nswap=self.para
		xlo,xhi,ylo,yhi,zlo,zhi,lx,ly,lz=box
		downP=upP;
		down11=xlo+deta*downP;
		down12=down11+deta;
		down22=xhi-deta*downP;
		down21=down22-deta;
		up12=xlo+lx/2-deta*upP;
		up11=up12-deta;
		up21=xlo+lx/2+deta*upP;
		up22=up21+deta;
		self.lp=up11-down11;
		print "region            up1    block   %f %f  INF INF INF INF  units box"%(up11,up12)
		print "region            up2    block   %f %f  INF INF INF INF  units box"%(up21,up22)
		print "region            up  union 2 up1 up2"
		print "region            down1  block  %f %f   INF INF INF INF units box"%(down11,down12)
		print "region            down2  block  %f %f   INF INF INF INF units box"%(down21,down22)
		print "region            down union 2 down1 down2"
	def renderCompute(self):
		print "compute           up_temp    all  temp/region up"
		print "compute           down_temp  all  temp/region down"
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
		print "fix getEqu  all  nvt temp %f %f %f"%(self.T,self.T,self.dtime)
	def renderElim(self):
		if self.nvt:
			print "fix getEqu  all  nvt temp %f %f %f"%(self.T,self.T,self.dtime)
		else:
			print "fix nve all nve"
	def renderFlux(self):
		print "fix	temp_profile    all    ave/spatial  1  %d  %d  x  lower  %f      v_temp  v_jx file  tempProfile.txt  norm sample units box"%(self.aveRate,self.aveRate,self.deta)
		# 输出热流空间分布,不计算热导率
		if self.jprofile==1:
			print "dump jprofile all custom %d jprofile.txt id v_jx v_jy v_jz v_temp v_jcx v_jcy v_jcz vx vy vz x y z"%(self.dumpRate)
			print "dump_modify  jprofile sort id"
	def renderSwap(self):
		lx=self.lx;
		lp=self.lp
		deta=self.deta
		Nswapbins=2*int(lx/(2*self.nswap*self.deta))
		factor=lx/(self.excRate*self.timestep)/(2*self.S)/(1/lp)*self.zfactor*self.tcfactor
		print "fix delta_out  all  ave/time  1  %d  %d  v_delta_temp   file  delta_temp.txt"%(self.aveRate,self.aveRate)
		print "fix heat_swap   all  thermal/conductivity  %d  x   %d"%(self.excRate,Nswapbins)
		print "fix e_exchange  all  ave/time  %d  %d  %d  f_heat_swap  file  swap.txt"%(self.excRate,self.excNum,self.aveRate)
		print "variable thermal_conductivity equal f_e_exchange/(1e-10+f_delta_out)*%f"%(factor)
		print "fix thermal_conductivity_out  all  ave/time  %d  1   %d  v_thermal_conductivity   file  thermal_conductivity.txt"%(self.aveRate,self.aveRate)
		
class ijDevice(mpDevice):
	def renderRegion(self):
		mpDevice.renderRegion(self)
		icoldl=self.xlo;
		icoldr=icoldl+self.nswap*self.deta;
		ihotl=self.xlo+self.lx/2;
		ihotr=ihotl+self.nswap*self.deta;
		print "region	hot block	%f %f   INF INF INF INF units box"%(ihotl,ihotr)
		print "region	cold block %f	%f	INF INF INF INF  units box"%(icoldl,icoldr)
	def renderSwap(self):
		factor=1/(2*self.S)/(1/self.lp)*self.zfactor*self.tcfactor
		swap=self.swapEnergyRate
		ave=self.aveRate
		exc=self.excRate
		print "fix delta_out  all  ave/time  1  %d  %d  v_delta_temp   file  delta_temp.txt"%(ave,ave)
		print "fix               hot   all  heat  %d   %f   region hot"%(exc,swap)
		print "fix               cold  all  heat  %d  -%f   region cold"%(exc,swap)
		print "variable          thermal_conductivity equal %f/(1e-10+f_delta_out)*%f"%(swap,factor)
		print "fix thermal_conductivity_out  all  ave/time  %d  1   %d  v_thermal_conductivity   file  thermal_conductivity.txt"%(ave,ave)
