from aces.materials.POSCAR import structure as Material
class structure(Material):
		def getPOSCAR(self):
			return """Ti2C
1.0
        3.0690870285         0.0000000000         0.0000000000
       -1.5345435143         2.6579073331         0.0000000000
        0.0000000000         0.0000000000        13.7363662720
    C     Ti
    1      2
Direct
     0.000000000         0.000000000         0.500000000
     0.666666985         0.333332986         0.583486021
     0.333332986         0.666666985         0.416514009

"""

"""
		def csetup(self):
				import numpy as np
				b=np.array([-0.25038,-0.968148,0.01,
						0.968148,-0.25038,0.01,
						0.01,0.01,1.0]).reshape([3,3])
				a=np.array([
						[0,1.0,0],
						[-1.0,-1.0,0],
						[0,0,1.0]
						])
				self.premitive=b.T.dot(a.T)
"""