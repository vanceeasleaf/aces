from aces.materials.POSCAR import structure as material
class structure(material):
	def getPOSCAR(self):
		return self.getMinimized()
		return """Mo  N                                   
	 1.0     
		 2.98   0    0
		 1.49   2.5807557   0
		 0    0   25
	 Mo   N 
		 1     2
Direct
	0.1666666666666643  0.6666666666666643  0.5000000000000000
	0.8333333333333357  0.3333333333333357   0.456
	0.8333333333333357  0.3333333333333357  0.544

"""
	def csetup(self):
			from ase.dft.kpoints import ibz_points
			#self.bandpoints=ibz_points['hexagonal']
			
			import numpy as np
			x=0.5*np.cos(np.arange(8)/8.0*2.0*np.pi)
			y=0.5*np.sin(np.arange(8)/8.0*2.0*np.pi)
			self.bandpath=['Gamma']
			for i in range(8):
				if(np.abs(x[i])>0.2):x[i]/=np.abs(x[i])*2.0
				if(np.abs(y[i])>0.2):y[i]/=np.abs(y[i])*2.0
				self.bandpoints['X'+str(i)]=[x[i],y[i],0.0]
				self.bandpath.append('X'+str(i))   
				self.bandpath.append('Gamma')    
	def getMinimized(self):
		return """POSCAR file written by OVITO
1.0
				2.9916000366         0.0000000000         0.0000000000
				1.4957998991         2.5908014232         0.0000000000
				0.0000000000         0.0000000000        25.0000000000
	 Mo    N
		1    2
Direct
		 0.000000000         0.000000000         0.500000000
		 0.666666687         0.666666687         0.455509961
		 0.666666687         0.666666687         0.544490039		
"""