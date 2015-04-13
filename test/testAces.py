#encoding: utf-8
import unittest
from aces import Aces,ThermalConductivity
class testAces(unittest.TestCase):
	def setUp(self):
		pass
	def tearDown(self):
		pass

	def testCommit(self):
		aces=Aces(stratege='uniform',procs=(1,4,'q1.4'))
		for i in xrange(31):
			units=Units('metal')
			opt={ method:'nvt'
				, units:units
				, species:'nvt-small'
				, thick:units.L(1.44)
				, langevin:0
				, usinglat:1
				, timestep:units.t.metal(.5e-3)
				, latx:11
				, runTime:10000000
				}
			app=ThermalConductivityApp(opt)
			aces.commit(app)
		aces.run() #产生所需的文件，记录所有参数到app.json
	
	def testQuery(self):
		aces=Aces(rebuild=True)
		opt={ upP:12}
		aces.query(opt)	
if __name__ =='__main__':
	unittest.main()
