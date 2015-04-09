#encoding: utf-8
import sys,os
print __file__
strfilepath=os.path.dirname(os.path.realpath(__file__))+"/../.."
sys.path.append(strfilepath)
import unittest

from src.profile import profile
class testProfile(unittest.TestCase):
	def setUp(self):
		pass
	def tearDown(self):
		pass

	def testGetTempProfile(self):
		profile().getTempProfile(1,12,2,5,5,5)
if __name__ =='__main__':
	unittest.main()
