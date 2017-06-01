from aces.materials.POSCAR import structure as Material
class structure(Material):
  def getPOSCAR(self):
    return """  findsym-output
1.0
        4.3336000443         0.0000000000         0.0000000000
       -2.1668000221         3.7530077282         0.0000000000
        0.0000000000         0.0000000000        25.3406600952
   Sn   Te
    2    2
Direct
     0.333333343         0.666666687         0.455769986
     0.666666627         0.333333313         0.544229984
     0.000000000         0.000000000         0.606950045
     0.000000000         0.000000000         0.393049985
"""
  def csetup(self):
    from ase.dft.kpoints import ibz_points
    self.bandpoints={'Gamma':[0,0,0],'Y':  [0,0.5,0],'M':[0.33, 0.33, 0],
        'L':[0.67, 0.33, 0],'X':[0.5, 0, 0 ]  }
    self.bandpath=['Gamma','Y','M','Gamma','L','X','Gamma']