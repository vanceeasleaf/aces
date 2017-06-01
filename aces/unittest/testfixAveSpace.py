#encoding: utf-8
import unittest
from aces.tools import write
from aces.tools import shell_exec
from aces.fixAveSpace import fixAveSpace
import numpy as np
class testfixAveSpace(unittest.TestCase):
	def setUp(self):
		content="""# Spatial-averaged data for fix temp_profile and group main
# Timestep Number-of-bins
# Bin Coord Ncount v_temp v_jx
100000 4
  1 14.961 0 0 0
  2 17.961 0 0 0
  3 20.961 10 291.241 -0.00395033
  4 23.961 10.9971 291.11 -0.00836294
200000 4
  1 14.961 0 0 0
  2 17.961 0 0 0
  3 20.961 10 287.894 0.00463008
  4 23.961 10.9962 291.343 0.00275501
300000 4
  1 14.961 0 0 0
  2 17.961 0 0 0
  3 20.961 10 290.52 0.0042
  4 23.961 10.9961 304.537 0.00489558
"""
		write(content,'tmpfixAveSpace.txt')
	def tearDown(self):
		shell_exec("rm tmpfixAveSpace.txt")
	def testinit(self):
		fas=fixAveSpace('tmpfixAveSpace.txt')
		self.assertEqual("Spatial-averaged data for fix temp_profile and group main",fas.title)
		self.assertEqual(fas.quants,["v_temp","v_jx"])
		self.assertEqual([fas.snapStep,fas.nbin,fas.nstep],[100000,4,3])
	def testgetIStep(self):
		fas=fixAveSpace('tmpfixAveSpace.txt')
		print fas.getIStep(1)[2][:,1]
		x=np.abs(fas.getIStep(1)[2][:,1]-np.array([0,0,0.00463008,0.00275501]))<0.01
		self.assertEqual(x.all(),True)
		x=np.abs(fas.getIStep(1)[2][:,0]-np.array([0,0,287.894,291.343]))<0.01
		self.assertEqual(x.all(),True)
		x=np.abs(fas.getIStep(1)[1]-np.array([0,0,10,10.9962]))<0.01
		self.assertEqual(x.all(),True)
	def testgetIbin(self):
		fas=fixAveSpace('tmpfixAveSpace.txt')
		x=np.abs(fas.getIbin(1)[:,1]-np.array([0,0,0]))<0.01
		self.assertEqual(x.all(),True)
		x=np.abs(fas.getIbin(3)[:,0]-np.array([291.11,291.343,304.537]))<0.01
		self.assertEqual(x.all(),True)
def test():
	unittest.main()
if __name__ =='__main__':
	test()
