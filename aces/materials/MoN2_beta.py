from aces.materials.POSCAR import structure as Material
class structure(Material):
	def getPOSCAR(self):
		return self.getMinimized()
		return """Mo  N                                   
	 1.0     
		 3.30   0    0
		 1.65   2.85788   0
		 0    0   25
	 Mo   N 
		 1     2
Direct
	0.1666666666666643  0.6666666666666643  0.5000000000000000
	0.8333333333333357  0.3333333333333357  0.4692
	0.8333333333333357  0.3333333333333357  0.5308

"""
	def csetup(self):
		from ase.dft.kpoints import ibz_points
		self.bandpoints=ibz_points['hexagonal']
		self.bandpath=['Gamma','M','K','Gamma']   			
	def getMinimized(self):
		return """POSCAR file written by OVITO
1.0
     3.3250620764572432    0.0000000000000000    0.0000000000000000
     1.6625310382286216    2.8795843074125806    0.0000000000000000
     0.0000000000000000    0.0000000000000000   25.0000000000000000
   Mo    N
    1    2
Direct
     0.000000000         0.000000000         0.500000000
     0.666666687         0.666666687         0.455509961
     0.666666687         0.666666687         0.544490039		
"""