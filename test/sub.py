#encoding: utf-8
from aces import Aces,Runner,NormalStrategy
from app.thermalConductivityApp import ThermalConductivityApp
from units import Units
Runner.run(sub())
class sub(Aces):
	def submit(self):
		self.strategy=NormalStrategy(procs=(1,4,'q1.4'))
		for i in xrange(31):
			units=Units('metal')
			opt={ method:'nvt'
				, units:units
				, species:'nvt-small'
				, thick:units.L(1.44)
				, langevin:0
				, usinglat:1
				, timestep:units.t.lj(.5e-3)
				, latx:11
				, runTime:10000000
				}
			app=ThermalConductivityApp(opt)
			self.register(app)

	def query(self):
		self.queryOption={ upP:12}
		self.want=['lax','kappa','ratio']

