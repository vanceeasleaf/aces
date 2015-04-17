#encoding:utf8
class gkDevice:
	def __init__(self,box,corRate,timestep,kb,T,zfactor,tcfactor,fourierTc,computeTc,gstart,corNum,aveRate,jcf,dtime):
		self.para=[box,corRate,timestep,kb,T,zfactor,tcfactor,fourierTc,computeTc,gstart,corNum,aveRate,jcf]
		self.T=T;self.dtime=dtime
	def renderRegion(self):
		pass
	def renderCompute(self):
		pass	
	def renderVariable(self):
		pass
	def renderEqu(self):
		print "fix getEqu  all  nvt temp %f %f %f"%(self.T,self.T,self.dtime)
	def renderElim(self):
		print "fix getEqu  all  nvt temp %f %f %f"%(self.T,self.T,self.dtime)
	def renderTemp(self):
		pass
	def renderFlux(self):
		box,corRate,timestep,kb,T,zfactor,tcfactor,fourierTc,computeTc,gstart,corNum,aveRate,jcf=self.para
		xlo,xhi,ylo,yhi,zlo,zhi,lx,ly,lz=box
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
	def renderSwap(self):
		pass