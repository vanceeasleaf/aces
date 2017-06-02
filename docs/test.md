# Test of the algorithm of apply TMD, strain

## knottmd

``` python
from aces import Aces
class sub(Aces):
	def submit(self):
		opt=dict(
			units="metal",
			species="graphene_knot",
			method="nvt",
			nodes=1,
			procs=12,
			queue="q1.1",
			runTime=10000000
			,runner="mdTc"
		)
		app=dict(latx=70,laty=2)
		self.commit(opt,app);
if __name__=='__main__':
	sub().run()
```
## largedt

``` python
from aces import Aces
class sub(Aces):
	def submit(self):
		opt=dict(
			units="metal",
			species="graphene_knot",
			method="nvt",
			nodes=1,
			procs=1,
			queue="q1.1",
			runTime=10000000
			,runner="mdTc",dT=100,useMini=False,laty=1,latz=1,latx=3
		)
		for i in range(3,20,4):
			app=dict(atomfile='POSCAR',timestep=.3e-3/i)
			self.commit(opt,app);
if __name__=='__main__':
	sub().run()
```
## testKnot
> test if aces could give correct knot

``` python
from aces import Aces
class sub(Aces):
	def submit(self):
		opt=dict(
			units="metal",
			species="graphene_knot",
			method="nvt",
			nodes=1,
			procs=12,
			queue="q1.1",
			runTime=10000000
			,runner="mdTc"
		)
		app=dict(latx=70,laty=2)
		self.commit(opt,app);
if __name__=='__main__':
	sub().run()
```

## strain-knot

``` python
from aces import Aces
class sub(Aces):
	def submit(self):
		opt=dict(
			units="metal",
			species="graphene_knot",
			method="nvt",
			nodes=1,
			procs=4,
			queue="q1.1",
			runTime=500000
			,runner="strain"
		)
		for T in range(100,300,20):
			app=dict(equTime=200000,T=T,strainStep=1000,
			maxStrain=0.3,timestep=.3e-3,latx=70,laty=2)
			self.commit(opt,app);
			app=app.copy()
			app['maxStrain']=-0.4
			self.commit(opt,app);
if __name__=='__main__':
	sub().run()
```
## strainknotV
> vStrain=true. using end v but not average distortion

``` python
from aces import Aces
class sub(Aces):
	def submit(self):
		opt=dict(
			units="metal",
			species="graphene_knot",
			method="nvt",
			nodes=1,
			procs=4,
			queue="q1.1",
			runTime=500000
			,runner="strain"
		)
		for T in range(100,300,20):
			app=dict(vStrain=True,equTime=200000,T=T,
			strainStep=1000,maxStrain=0.3,timestep=.3e-3,latx=70,laty=2)
			self.commit(opt,app);
			app=app.copy()
			app['maxStrain']=-0.4
			self.commit(opt,app);
if __name__=='__main__':
	sub().run()
```
## testairebo
> text airebo potential, however knot does not exist

``` python
from aces import Aces
class sub(Aces):
	def submit(self):
		opt=dict(
			units="metal",
			species="graphene_knot",
			method="nvt",
			nodes=1,
			procs=12,
			queue="q1.1",
			runTime=10000000
			,runner="mdTc"
		)
		app=dict(airebo=True,timestep=.182e-3,latx=70,laty=2)
		self.commit(opt,app);
if __name__=='__main__':
	sub().run()
```
