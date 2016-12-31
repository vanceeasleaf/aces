from aces.materials.POSCAR import structure as material
class structure(material):
		def getPOSCAR(self):
			return """ACES POSCAR
1.0
3.61414867263 0.0 0.0
0.0 3.61414867263 0.0
0.0 0.0 13.1103467941
Co S H
2 2 2
Direct
0.5 0.0 0.52117306
0.0 0.5 0.52117306
0.5 0.5 0.614274831341
0.0 0.0 0.428071318659
0.5 0.5 0.714274831341
0.0 0.0 0.328071318659
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