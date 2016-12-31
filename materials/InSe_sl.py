from aces.materials.POSCAR import structure as material
class structure(material):
  def getPOSCAR(self):
		return """POSCAR file written by OVITO
1.0
        4.0803017616         0.0000000000         0.0000000000
       -2.0003423035         3.5563319623         0.0000000000
       -0.0888649683         0.0519736559        27.5121021594
   In   Se
    2    2
Direct
     0.331981868         0.668018103         0.558396757
     0.331082582         0.668917298         0.455985188
     0.995835781         0.004164212         0.409338266
     0.998198986         0.001800989         0.605088770
"""
  def csetup(self):
    from ase.dft.kpoints import ibz_points
    self.bandpoints=ibz_points['hexagonal']
    self.bandpath=['Gamma','K','M','Gamma']
			
		