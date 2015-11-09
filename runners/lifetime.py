#encoding:utf8
from aces.runners import Runner
from aces.runners.correlation import runner as Crun
from aces.runners.phonopy import runner as Prun
from aces.graph import plot,imshow,scatter
import numpy as np
from  aces.tools import *
hbar=6.6260755e-34/3.14159/2.0
kb=1.3806488e-23
#Bose-Einstein Distribution
def BE(w,T):
	w=np.array(w)
	t= hbar*w/kb/T
	#return np.exp(-t)
	return 1.0/(np.exp(t)-1.0)
class runner(Runner):
	def generate(self):
		prun=Prun(self.m)
		prun.run()
		self.corr()
		self.fc()
		self.sed()
		#self.band()
		#self.drawband()
	def corr(self):
		crun=Crun(self.m)
		crun.run()
		#crun.vd.life_yaml(correlation_supercell=self.m.correlation_supercell)
		#self.drawlifetime()
	def fc(self):
		c=self.m.correlation_supercell
		q=[]
		for i in range(c[0]):
			for j in range(c[1]):
				for k in range(c[2]):
					q.append([float(i)/c[0],float(j)/c[1],float(k)/c[2]])
		m=self.m
		
		mkcd('qpoints')
		cp('../FORCE_CONSTANTS','.')
		cp('../disp.yaml','.')
		cp('../POSCAR','.')
		Prun(m).getqpoints(q)
		cd('..')
	def dos(self):
		crun=Crun(self.m)
		crun.dos()
	def nmaq(self,k=[0,0,0],test=False):
		from aces.runners.vdos import vdos
		correlation_supercell=self.m.correlation_supercell
		vdos(self.m.timestep).lifenmaq(k=k,correlation_supercell=correlation_supercell,test=test)
	def nma(self):
		from aces.runners.vdos import vdos
		correlation_supercell=self.m.correlation_supercell
		vdos(self.m.timestep).lifenma(correlation_supercell=correlation_supercell)
	def sed(self):
		from aces.runners.vdos import vdos
		correlation_supercell=self.m.correlation_supercell
		vdos(self.m.timestep).lifesed(correlation_supercell=correlation_supercell)
	def band(self):
		from aces.runners.vdos import vdos
		correlation_supercell=self.m.correlation_supercell
		vdos(self.m.timestep).sed_band(correlation_supercell=correlation_supercell)
	def drawband(self):
		sed=np.load('sed.npy')
		n=len(sed)
		x=np.arange(n)
		w=np.load('wtick.npy')
		filter=w<60
		sed=sed[:,filter]
		#X,W=np.meshgrid(x,w)
		#scatter(X,W,sed,'Wave Vector','Frequency (THz)','sed.png',marker_size=2,marker='s')
		imshow(np.log(sed.T),'sed.png',extent=[0,1,0,1])
	def drawlifetime(self):
		a=np.loadtxt('life.txt',skiprows=1)
		om=a[:,3]
		tao=np.abs(1/a[:,5])
		n=len(tao)
		filter=tao<1e10
		om=om[filter]
		tao=tao[filter]
		plot((om,'Frequency (THz)'),(tao,'Relaxation Time (ps)'),'tao_freq.png',grid=True,scatter=True,logx=True,logy=True)
		
		
		v=np.loadtxt('groupv/v.txt')[:n,1]
		v=v[filter]
		l=v*tao
		to_txt(['freq','lamda'],np.c_[om[om.argsort()],l[om.argsort()]],'lamda_freq.txt')
	
		plot((om,'Frequency (THz)'),(l,'Mean Free Path (Angstrom)'),'lamda_freq.png',grid=True,scatter=True,logy=True)
		T=self.m.T
		w=om*1e12
		c=hbar*w*(BE(w,T+0.005)-BE(w,T-0.005))*100.0
		to_txt(['freq','capacity'],np.c_[om[om.argsort()],c[om.argsort()]],'capacity_freq.txt')
		plot((om[om.argsort()],'Frequency (THz)'),(c[om.argsort()],'Mode Specific Heat (J/K)'),'capacity_freq.png',grid=True,linewidth=2)
		to_txt(['freq','cumsumcapacity'],np.c_[om[om.argsort()],c[om.argsort()].cumsum()],'cumsumcapacity_freq.txt')
		plot((om[om.argsort()],'Frequency (THz)'),(c[om.argsort()].cumsum(),'Acummulate Specific Heat (J/K)'),'cumsumcapacity_freq.png',grid=True,linewidth=2)

		from ase import io 
		atoms=io.read('POSCAR')
		V=np.linalg.det(atoms.cell)
		k=l*c*v/V*1e12*1e10
		to_txt(['freq','kappa'],np.c_[om[om.argsort()],k[om.argsort()]],'kappa_freq.txt')
		plot((om,'Frequency (THz)'),(k,'Mode Themal Conductivity (W/mK)'),'kappa_freq.png',grid=True,scatter=True)

		to_txt(['freq','cumsumkappa'],np.c_[om[om.argsort()],k[om.argsort()].cumsum()],'cumsumkappa_freq.txt')
		plot((om[om.argsort()],'Frequency (THz)'),(k[om.argsort()].cumsum(),'Acummulate Themal Conductivity (W/mK)'),'cumsumkappa_freq.png',grid=True,linewidth=2)

		to_txt(['lamda','cumsumkappa'],np.c_[l[l.argsort()],k[l.argsort()].cumsum()],'cumsumkappa_lamda.txt')
		plot((l[l.argsort()],'Mean Free Path (Angstrom)'),(k[l.argsort()].cumsum(),'Acummulate Themal Conductivity (W/mK)'),'cumsumkappa_lamda.png',grid=True,linewidth=2)
