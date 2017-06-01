# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2016-09-06 16:41:54
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-01 21:40:29
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
		pass#['Gamma','Y','T','X','Gamma']

	def setup(self):
		self.forceThick=False
		self.elements=['Zr','Ni','Sn']
		self.premitive/=np.array([self.latx,self.laty,self.latz])

	def lmp_structure(self):
		pos=np.array([
			[0.5000000000000000,   0.5000000000000000,   0.500000000000000],
			[0.2500000000000000,   0.2500000000000000,   0.250000000000000],
			[0.0000000000000000,   0.0000000000000000,   0.000000000000000]
		])
		cell=6.11*np.array([
			[0.0, 0.5, 0.5],
			[0.5, 0.0, 0.5],
			[0.5, 0.5, 0.0]
		])
		atoms = Atoms('ZrNiSn',scaled_positions=pos, cell=cell)
		atoms=atoms.repeat([self.latx,self.laty,self.latz])
		atoms.set_pbc([self.xp,self.yp,self.zp])
		return atoms
			
		