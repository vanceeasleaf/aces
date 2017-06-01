from aces.materials  import Material
from aces.modify import get_unique_atoms
from ase import Atoms,Atom
from math import pi,sqrt,atan
import numpy as np
class structure(Material):
	def set_parameters(self):
		pass
	def setup(self):
		pass
	def lmp_structure(self):		
		col=self.unitcell(self.laty,self.latx)		
		col.set_pbc([self.xp,self.yp,self.zp])
		atoms=get_unique_atoms(col)
		cell=atoms.cell*self.bond
		atoms.set_cell(cell,scale_atoms=True)
		atoms.center()
		return atoms
		
	
	def unitcell(self,latx,laty):
		pos2=np.array([1.02325,3.54463,0
,0.34109,1.18154,0
,5.7984,4.72618,0
,5.11624,2.36309,0
,7.50382,2.95387,0
,3.41082,4.13541,0
,2.72866,1.77231,0
,8.18599,5.31696,0
,6.82166,0.59078,0
,2.04649,14.76931,0
,8.52706,14.17854,0
,5.45733,11.22468,0
,3.75191,12.997,0
,1.36433,12.40622,0
,0.68216,10.04313,0
,6.13949,13.58776,0
,3.06975,10.6339,0
,7.8449,11.81546,0
,6.48057,7.08926,0
,4.77515,8.86159,0
,7.16273,9.45235,0
,1.70542,5.90772,0
,4.09299,6.4985,0
,2.38758,8.27081,0
,0,7.68004,0
,4.43408,0,0
,7.16273,4.33233,0
,4.77515,3.74156,0
,2.38758,3.15078,0
,5.45733,6.10465,0
,3.06975,5.51387,0
,0.68216,4.92311,0
,7.8449,6.69542,0
,4.09299,1.37846,0
,1.70542,0.7877,0
,8.18599,0.19692,0
,6.48057,1.96924,0
,3.41082,14.37546,0
,1.02325,13.7847,0
,0.34109,11.42161,0
,2.72866,12.01237,0
,5.7984,14.96624,0
,5.11624,12.60315,0
,7.50382,13.19392,0
,8.52706,9.0585,0
,6.82166,10.83083,0
,3.75191,7.87696,0
,6.13949,8.46774,0
,1.36433,7.2862,0
,2.04649,9.64928,0
,0,2.56002,0
,4.43408,10.24006,0]).reshape(-1,3)/1.42
		pos1=[6.928400,13.000369,0.000000
			,7.794450,16.500469,0.000000]
		phi=pi/2-atan((pos1[4]-pos1[1])/(pos1[3]-pos1[0]))
		cbond=np.linalg.norm((pos1[4]-pos1[1],pos1[3]-pos1[0],0))
		dx=sqrt(3)*cbond;
		dy=3*cbond;
		atoms=Atoms()

		for i,coord in enumerate(pos2):
			ele=['C','N'][i<26]
			atom=Atom(ele,coord)
			atoms.append(atom)
		#atoms.rotate('z',phi)
		atoms.set_cell([dx,dy,10.0])
		atoms=atoms.repeat((2,2,1))
		col=atoms.repeat((latx,laty,1))
		return col
		
			
		