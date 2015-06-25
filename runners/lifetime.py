#encoding:utf8
from aces.runners import Runner
from aces.runners.correlation import runner as Crun
from aces.runners.phonopy import runner as Prun
class runner(Runner):
	def generate(self):
		crun=Crun(self.m)
		prun=Prun(self.m)
		prun.run()
		crun.run()
		crun.vd.life_yaml(correlation_supercell=self.m.correlation_supercell)