#encoding: utf-8
import unittest
from aces.tools import write,shell_exec
from aces.pizza.data import data
from aces.lammpsdata import lammpsdata
from math import pi,sqrt,cos
from ase import Atoms

import numpy as np
from numpy import allclose as cc
class testlammpsdata(unittest.TestCase):
	def setUp(self):
		pass

	def tearDown(self):
		pass
	def testmergeVec(self):
		d=lammpsdata(Atoms())
		direct,phi=d.mergeVec([1,0,0],[0,2,0])
		self.assertTrue(cc(direct,[0,0,1]))
		self.assertEqual(phi,pi/2)
		direct,phi=d.mergeVec([1,0,1],[1,1,1])
		self.assertTrue(cc(direct,[-1/sqrt(2),0,1/sqrt(2)]))
		self.assertEqual(cos(phi),sqrt(2.0/3))
		
		direct,phi=d.mergeVec([1,0,0],[1,0,0])
		self.assertTrue(cc(direct,[1,0,0]))
		self.assertEqual(cos(phi),1)
	def testgetTypes(self):
		a=Atoms('CNC',[[0,0,1]]*3)
		d=lammpsdata(a)
		self.assertEqual(['C','N'],d.getTypes())
		a=Atoms('CNSiCCSiSBi',[[0,0,1]]*8)
		d=lammpsdata(a)
		self.assertEqual(['C','N','Si','S','Bi'].sort(),d.getTypes().sort())
	def testget_rotated_atoms(self):
		a=Atoms('C2',[[1,0,1],[0,1,1]])
		a.set_cell([[1,0,1],[0,1,1],[0,0,1]])
		d=lammpsdata(a)
		unit=d.get_rotated_atoms()
		self.assertTrue(cc(unit.cell[:2],[[sqrt(2),0,0],[sqrt(2)/2,sqrt(2)*sqrt(3)/2,0]]))

	def testwritedata(self):
		a=Atoms('C2',[[1,0,1],[0,1,1]])
		a.set_cell([[1,0,1],[0,1,1],[0,0,1]])
		d=lammpsdata(a,['C','N'])
		d.writedata()
		b=data('structure')
		self.assertEqual(b.headers["atoms"],2)
		self.assertEqual(b.headers["atom types"],2)
		self.assertEqual(b.title ,'C2\n')
		self.assertTrue(cc(b.headers["xlo xhi"],[0,sqrt(2)]))
		self.assertTrue(cc(b.headers["ylo yhi"],[0,sqrt(2)*sqrt(3)/2]))
		self.assertTrue(cc(b.headers["xy xz yz"][0:1],[sqrt(2)/2]))
		b.map(1,'id',2,'type',3,'x',4,'y',5,'z')
		atoms=b.viz(0)[2]
		self.assertTrue(cc(atoms,[[1,1,sqrt(2),0,0],[2,1,sqrt(2)/2,sqrt(2)*sqrt(3)/2,0]]))
		shell_exec('rm structure')
def test():
	unittest.main()
if __name__ =='__main__':
	test()
