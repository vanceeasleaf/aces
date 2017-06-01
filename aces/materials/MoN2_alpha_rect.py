from aces.materials.POSCAR import structure as Material
class structure(Material):
	def getPOSCAR(self):
		return self.getMinimized()
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
			#self.bandpath=['Gamma',"X2"]
	def getMinimized(self):
		return """Mo  N 
 1.0000000000000000
     2.9916000366000000    0.0000000000000000    0.0000000000000000
     0.0000000000000000    5.1814560994168932    0.0000000000000000
     0.0000000000000000    0.0000000000000000   25.0000000000000000
  Mo   N
   2   4
Direct
  0.5000000000000000  0.5000000000000000  0.5000000000000000
  0.0000000000000000  0.0000000000000000  0.5000000000000000
  0.5000000000000000  0.8333333333333335  0.4555099610000000
  0.5000000000000000  0.8333333333333335  0.5444900390000000
  0.0000000000000000  0.3333333333333333  0.4555099610000000
  0.0000000000000000  0.3333333333333333  0.5444900390000000

"""