from aces.materials  import Material
from aces.modify import get_unique_atoms
from ase import Atoms,Atom
from math import pi,sqrt
from ase.dft.kpoints import ibz_points
from aces import config
from ase.lattice import bulk
import numpy as np
from aces.tools import *
class structure(Material):
	def set_parameters(self):
		self.modulation=False
		self.cu=True
		self.se=False;
		pass#['Gamma','Y','T','X','Gamma']

	def setup(self):
		self.forceThick=False
		Q='Te'
		if self.se:
			Q='Se'
		self.elements=['Sn',Q]
		#self.bandpoints=ibz_points['fcc']
		#self.bandpath=['Gamma','K','X','Gamma','L']
		self.premitive/=np.array([self.latx,self.laty,self.latz])

	def lmp_structure(self):
		Q='Te'
		if self.se:
			Q='Se'
		pos=np.array([[0,0,0],[.5,.5,.5]])
		cell=3.2019243179955601*2.0*np.array([[0,.5,.5],[.5,0,.5],[.5,.5,0]])
		if self.modulation:
			pos=np.array([[0.0068692286303052 , 0.0063691994833112,  0.0050618520044743],
						  [0.4836396400495982 , 0.5068917001717190 , 0.5092467886015576]])
			cell=np.array([
				[-0.0028854017661314  ,  3.2115195044645874 ,   3.2127505583949509],
			    [ 3.2011692839678605  ,  0.0049602662391963 ,   3.2049160494026143],
			    [ 3.2013927516914373  ,  3.2038936635259327 ,   0.0047286386862793]])

		atoms = Atoms('Sn'+Q,scaled_positions=pos, cell=cell)
		if(self.cu):
			pos=np.array([
						[0.0104073240749400,  0.0121251794234709,  0.0076458899654810],
						[0.0104073240749400,  0.5121251794234709,  0.5076458899654810],
						[0.5104073240749400,  0.0121251794234709,  0.5076458899654810],
						[0.5104073240749400,  0.5121251794234709,  0.0076458899654810],
						[0.4997298381725912,  0.5014483020557208,  0.4969676463570991],
						[0.4997298381725912,  0.0014483020557208,  0.9969676463570991],
						[0.9997298381725912,  0.5014483020557208,  0.9969676463570991],
						[0.9997298381725912,  0.0014483020557208,  0.4969676463570991]
				])
			cell=np.array([
				[6.4134466126604810,    0.0123852356641863,    0.0130369683911169],
     			[0.0132818917551842,    6.4134441681750927,    0.0132871756523024],
     			[0.0126358349919128,    0.0123870727843269,    6.4134474786386608]
     			])
			atoms = Atoms('Sn4'+Q+'4',scaled_positions=pos, cell=cell)
		atoms=atoms.repeat([self.latx,self.laty,self.latz])
		atoms.set_pbc([self.xp,self.yp,self.zp])
		return atoms
			
		