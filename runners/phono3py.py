#encoding:utf8
from aces.tools import *
import aces.config as config
from ase.io import read
from ase.io.vasp import write_vasp
from aces.binary import pr
from aces.runners import Runner
from aces.graph import plot,series
import numpy as np
from aces.runners.phonopy import runner as Runner
class runner(Runner):
			
	def force_constant(self,files):
		cmd=config.phono3py+"--cf3 "
		for file in files:
			dir="dirs/dir_"+file
			cmd+=dir+'/vasprun.xml '
		write(cmd,'post.sh')
		#generate FORCE_SETS
		passthru("sh post.sh")
		m=self.m
		print "Create fc2.hdf and fc3.hdf"
		passthru(config.phono3py+"-c POSCAR --dim='%s'"%(m.dim))
		
	def generate_supercells(self):
		m=self.m
		#generate supercells
		passthru(config.phono3py+"-d --dim='%s' -c POSCAR"%(m.dim))

		
	def generate(self):
		m=self.m
		self.minimizePOSCAR()
		self.generate_supercells()
		files=shell_exec("ls *-*").split('\n')
		self.getvasprun(files)
		self.force_constant(files)
		if m.pho3bte:
			self.sheng()
		else:
			self.caltc()
		#Thermal conductivity calculation
	def caltc(self):
		m=self.m
		passthru(config.phono3py+'  --dim="'+m.dim+'" -c POSCAR --mesh="'+' '.join(map(str,m.kpoints))+'" --fc3 --fc2  --thm --br')
		filename="kappa-m%s.hdf5"%''.join(map(str,m.kpoints))
		passthru(config.pypath+'kaccum   --mesh="'+' '.join(map(str,m.kpoints))+'"  POSCAR '+filename+' |tee kaccum.dat')
	def sheng(self):
		
		self.writeFC()
		m=self.m
		from  aces.runners.shengbte import runner as shengbte
		a=shengbte(m)
		mkcd('SHENG')
		cp('../FORCE_CONSTANTS_2ND','.')
		cp('../FORCE_CONSTANTS_3RD','.')
		a.getControl()
		#Thermal conductivity calculation
		print "START SHENGBTE..."
		#passthru(config.mpirun+" %s "%m.cores+config.shengbte)
		passthru(config.shengbte)
	def matrixFormat(self,mat):
		n,m=mat.shape
		s=""
		for k in range(n):
			for p in range(m):
				s+="\t%f"%mat[k,p]
			s+="\n"	
		return s
	def matrix3Format(self,mat):
		n,m,q=mat.shape
		s=""
		for k in range(n):
			for p in range(m):
				for t in range(q):
					s+="%d %d %d %f\n"%(k+1,p+1,t+1,mat[k,p,t])
		return s		
	def writeFC(self):
		print "writing text FORCE CONSTANTS 2 from hdf5"
		import h5py
		f=h5py.File('fc2.hdf5')
		fc2=f['fc2']
		natom=len(fc2)
		s="%d\n"%natom
		for i in range(natom):
			for j in range(natom):
				s+="%d\t%d\n"%(i+1,j+1)
				s+=self.matrixFormat(fc2[i,j])
		write(s,'FORCE_CONSTANTS_2ND')
		f.close()
		print "writing text FORCE CONSTANTS 3 from hdf5"
		f=h5py.File('fc3.hdf5')
		#phono3py phi is bac , sheng bte phi is bca and shengbte thirdorder need 100.0 scale
		fc3=np.einsum(f['fc3'],[0,2,1,3,5,4])/100.0
		natom=len(fc3)
		s=""
		n=0
		for i in range(natom):
			for j in range(natom):
				for k in range(natom):
					if np.allclose(np.zeros([3,3,3]),fc3[i,j,k],atol=1e-04):continue
					n+=1
					s+="\n%d\n"%n
					s+=self.maps2p(i,j,k)
					s+=self.matrix3Format(fc3[i,j,k])
					
		write("%d\n"%n+s,'FORCE_CONSTANTS_3RD')	
	def getS2p(self):
		atoms=read('POSCAR')
		import yaml
		data = yaml.load(open('disp_fc3.yaml').read())
		cell=np.array(data['lattice'])
		cellp=atoms.cell
		satoms=data['atoms']
		#symbols=''.join([o['symbol'] for o in satoms]])
		posp=atoms.get_scaled_positions()
		pos=np.array([o['position'] for o in satoms])
		#superatoms=Atoms(symbols,pos,scaled_positions=pos,cell=cell)
		#which cell is i j in
		vpos=pos*np.linalg.norm(cell,axis=1)/np.linalg.norm(cellp,axis=1)
		v=np.floor(vpos)
		
		celloffset=v.dot(cellp)
		vpos-=v
		s2p=-np.ones(len(vpos),dtype='int')
		for pd,p in enumerate(vpos):
			for id,a in enumerate(posp):
				if np.allclose(p,a):
					s2p[pd]=id+1
					break

		return s2p,celloffset
	def maps2p(self,i,j,k):
		if not hasattr(self, 's2p'):
			self.s2p,self.celloffset=self.getS2p()
		m=self.m
		v=self.celloffset

		s=m.toString(v[j]-v[i])+'\n'
		s+=m.toString(v[k]-v[i])+'\n'
		s+='%d %d %d\n'%(self.s2p[i],self.s2p[j],self.s2p[k])
		return s




