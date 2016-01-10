#encoding:utf8
class gkDevice:
	def __init__(self,hook,m):
		self.__dict__=dict(self.__dict__,**m.__dict__)
		hook.addAction('equ',self.renderEqu)
		hook.addAction('elimination',self.renderElim)
		hook.addAction('flux',self.renderFlux)
	def renderEqu(self):

		print "fix getEqu  all  nvt temp %f %f %f"%(self.T,self.T,self.dtime)
	def renderElim(self):
		if self.nvt:
			print "fix getEqu  all  nvt temp %f %f %f"%(self.T,self.T,self.dtime)
		else:
			print "fix nve all nve"
		print "fix mom all momentum 100 linear 1 0 0"
	def renderFlux(self):
		xlo,xhi,ylo,yhi,zlo,zhi,lx,ly,lz=self.box
		v=lx*ly*lz
		factor=self.corRate*self.timestep/(v*self.kb*self.T*self.T)*self.zfactor*self.tcfactor
		#/* 用傅里叶变换热流来计算热流关联函数*/
		if self.fourierTc:
			print "fix j_out all ave/time 1 1 1 c_jflux[1] c_jflux[2] c_jflux[3]  file  jin.txt "
			#print "dump vf all  custom 1 dump.vf id vx vy vz fx fy fz  "
			#print "dump_modify vf sort id"
		elif self.computeTc:
		#/* 用分子模拟论坛的compute扩展来计算热流关联函数,比lammps自带的更精确*/
			rfactor=self.tcfactor*self.zfactor
			print "variable          factor_ac equal 1.0"
			print "variable          factor_tc equal %f"%(rfactor)
			print "compute           tc all tc c_thermo_temp c_jflux v_factor_ac v_factor_tc x first %d 0 500000 %d"%(self.aveRate/self.corRate,self.corRate)
			print "fix               tc_out  all  ave/time  1  1  %d  c_tc   file  kappa.txt"%(self.corRate)
		else:
			corNum=self.aveRate/self.corRate
			print "fix ss all ave/correlate %d  %d  %d c_jflux[1] c_jflux[2] c_jflux[3] type auto ave running"%(self.corRate,corNum,self.aveRate)
			print "variable k11 equal trap(f_ss[3])*%f"%(factor)
			print "variable k22 equal trap(f_ss[4])*%f"%(factor)
			print "variable k33 equal trap(f_ss[5])*%f"%(factor)
			print "fix output all ave/time 1  1 %d v_k11  v_k22  v_k33 file kappa.txt"%(self.aveRate)
			#/* 定时输出热流自关联函数*/
			if self.jcf:
				print "fix out all ave/time %d  1   %d f_ss[3] f_ss[4] f_ss[5] mode vector file jcf.txt"%(self.aveRate,self.aveRate)
