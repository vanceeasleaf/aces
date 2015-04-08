#encoding: utf-8
import sys
sys.path.append("..")
from profile import profile
class testProfile(unittest.TestCase):
	def setUp(self):
		pass
	def tearDown(self):
		pass

	def testGetTempProfile(self):
		profile().getTempProfile(20,12,2,5,5,5)
if __name__ =='__main__':
	unittest.main()
