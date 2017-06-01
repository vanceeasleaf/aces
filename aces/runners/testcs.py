#encoding: utf-8
import unittest
from aces.runners.cs import runner
import numpy as np
class case1(unittest.TestCase):
	def setUp(self):
		# a trainsets of 2 set of 2 atoms 
		F=np.array(
			[[[1.0,2.0,3.0],
		    [4.0,5.0,6.0]],
		   [[-1.0,-2.0,-3.0],
		    [-4.0,-5.0,-6.0]]])
		u=np.array(
			[[[.3,.2,.1],
		    [.4,.3,.2]],
		   [[-2.0,-1.0,-0.5],
		    [-0.3,0.5,0.9]]])
		a=runner(NAH=3)
		a.getTrainSets(F,u)
		self.a=a
		self.F,self.u=F,u
	def tearDown(self):
		pass
	def testforceFile(self):
		pass
	def testgetTrainSets(self):
	
		a=self.a
		self.assertEqual(a.L,2)
		self.assertEqual(a.natom,2)
		self.assertEqual(a.rowr,[0,1,7,43])
	def testgetMatrix(self):
		a=self.a
		F,A=a.getMatrix(self.F,self.u)
		self.assertTrue(np.allclose(F,[[1.0,2.0,3.0,4.0,5.0,6.0],
		   [-1.0,-2.0,-3.0,-4.0,-5.0,-6.0]]))
		u1=[.3,.2,.1,.4,.3,.2]
		v1=np.array(u1)
		u2=[-2.0,-1.0,-0.5,-0.3,0.5,0.9]
		v2=np.array(u2)
		v3=[]
		for i in u1:
			for j in u1:
				v3.append(i*j)
		v3=np.array(v3)
		v4=[]
		for i in u2:
			for j in u2:
				v4.append(i*j)
		v4=np.array(v4)
		self.assertTrue(np.allclose(A,[np.hstack([-1,-v1,-v3/2]),np.hstack([-1,-v2,-v4/2])]))
	def testmulU(self):
		a=self.a
		x=np.array([1.0,2,3])
		self.assertTrue(np.allclose(a.mulU(x,1),[1.0,2,3]))
		self.assertTrue(np.allclose(a.mulU(x,2),np.array([1.0,2,3,2,4,6,3,6,9])/2))
	def testrebuild(self):
		a=self.a
		B=np.arange(43*6).reshape([43,6])
		phi=a.rebuild(B)
		self.assertTrue(np.allclose(phi[0],[[0.0,1,2],[3,4,5]]))
		self.assertTrue(np.allclose(phi[1][0][0],[[6,12,18],[7,13,19],[8,14,20]]))
		self.assertTrue(np.allclose(phi[1][1][0],[[9,15,21],[10,16,22],[11,17,23]]))
def test():
	unittest.main()
if __name__ =='__main__':
	test()
