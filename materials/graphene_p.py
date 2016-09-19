from aces.materials.POSCAR import structure as material
class structure(material):
	def getPOSCAR(self):
		return """ACES POSCAR                             
   1.00000000000000     
     2.4662118682569303   -0.0000000000000108    0.0000000000000002
     1.2331059341234423    2.1358021290232667   -0.0000000000000001
    -0.0000000000000007    0.0000000000000007    7.9507129814977562
   C 
     2
Direct
  0.1666666666669983  0.6666666666669983  0.5000000000000000
  0.8333333333330017  0.3333333333330017  0.5000000000000000

"""
	def csetup(self):
		from ase.dft.kpoints import ibz_points
		self.bandpoints=ibz_points['hexagonal']
		self.bandpath=['Gamma','M','K','Gamma']		
		