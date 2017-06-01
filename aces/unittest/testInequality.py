#encoding: utf-8
import unittest
import inequality as ie
class testInequality(unittest.TestCase):
	def setUp(self):
		self.ie=ie.inequality()

	def tearDown(self):
		self.ie=None

	def testreaddump(self):
		self.assertEqual(myclass.sum(1, 2), 3, 'test sum fail')
	def testsub(self):
		self.assertEqual(myclass.sub(2, 1), 1, 'test sub fail')   
		
if __name__ =='__main__':
	unittest.main()
