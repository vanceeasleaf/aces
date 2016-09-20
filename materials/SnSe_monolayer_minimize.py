from aces.materials.POSCAR import structure as material
class structure(material):
  def getPOSCAR(self):
		return """ACES POSCAR                             
   1.00000000000000     
     4.3715553598745540    0.0000000000000000    0.0000000000000000
     0.0000000000000000    4.2918241965252291    0.0000000000000000
     0.0000000000000000    0.0000000000000000   22.3387621007543906
   Sn   Se
     2     2
Direct
  0.0214617708073111  0.2500000000000000  0.4384455815054507
  0.5214617708069064  0.7500000000000000  0.5615544184945497
  0.4785382291930940  0.7500000000000000  0.4399051558263891
  0.9785382291930936  0.2500000000000000  0.5600948441736112
"""

  def csetup(self):
    from ase.dft.kpoints import ibz_points
    self.bandpoints=ibz_points['orthorhombic']
    self.bandpoints['T']=self.bandpoints['S']
    self.bandpath=['Gamma','Y','T','X','Gamma']
			
		