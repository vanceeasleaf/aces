#encoding: utf-8
import unittest
from aces.tools import write
from aces.tools import shell_exec
from aces.lineManager import lineManager
class testlineManager(unittest.TestCase):
	def setUp(self):
		content="""# Spatial-averaged data for fix temp_profile and group main
# Timestep Number-of-bins
# Bin Coord Ncount v_temp v_jx
100000 19
  1 14.961 0 0 0
  2 17.961 0 0 0
  3 20.961 10 291.241 -0.00395033
"""
		write(content,'tmptestlinemanager.txt')
	def tearDown(self):
		shell_exec("rm tmptestlinemanager.txt")
	def testparse(self):
		lm=lineManager('tmptestlinemanager.txt')
		self.assertEqual(7,lm.nline)
	def testgetLine(self):
		lm=lineManager('tmptestlinemanager.txt')
		self.assertEqual("# Bin Coord Ncount v_temp v_jx",lm.getLine(2))
		self.assertEqual("2 17.961 0 0 0",lm.getLine(5))

		
def test():
	unittest.main()
if __name__ =='__main__':
	test()
