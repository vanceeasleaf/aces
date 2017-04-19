import unittest
import myclass
class mytest(unittest.TestCase):
    def setUp(self):
        pass
    def tearDown(self):
        pass
    def testsum(self):
        self.assertEqual(myclass.sum(1, 2), 3, 'test sum fail')
    def testsub(self):
        self.assertEqual(myclass.sub(2, 1), 1, 'test sub fail')   
if __name__ =='__main__':
    unittest.main()
