#encoding: utf-8
import unittest
from ase import Atoms
from ase.io import read
from ase.calculators.neighborlist import NeighborList
from numpy import dot
from math import atan2,pi
class testAse(unittest.TestCase):
	def setUp(self):
		pass
	def tearDown(self):
		pass

	def testRead(self):
		atoms=read('range',format='lammps')
		atoms. set_pbc((1, 1, 1))
		sbs=atoms.get_chemical_symbols()
		c=[]
		for sb in sbs:
			if sb=='H':
				c.append('C')
			else:
				c.append('N')
				
		atoms.set_chemical_symbols(c)
	#	print atoms.get_chemical_symbols()	
		print atoms;
	def testNeighborlist(self):
		atoms=read('range',format='lammps')
		atoms. set_pbc((1, 1, 1))
		nl = NeighborList([0.8 for atom in atoms],self_interaction=False,bothways=True)
		nl.update(atoms)
		ang=[]		
		for i in xrange(3):
			indices, offsets = nl. get_neighbors(i)
			angs=[]
			for j, offset in zip(indices, offsets):
				pos= atoms. positions[j] + dot(offset, atoms. get_cell())-atoms.positions[i]
				ang1=atan2(pos[1],pos[0])+pi
				angs.append((j,ang1))
			newangs=sorted(angs,key=lambda d:d[1])
			print newangs
if __name__ =='__main__':
	unittest.main()
