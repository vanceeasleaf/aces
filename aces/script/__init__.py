# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2015-10-14 20:43:32
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-06 16:49:42
from aces.tools import *
import numpy as np
from ase import io
from aces.graph import pl,fig
from numpy.linalg import norm 
import itertools

def parseVasprun(vasprun,tag="forces"):
	forces = []
	for event, element in vasprun:
		if element.attrib['name'] == tag:
			for v in element.xpath('./v'):
				forces.append([float(x) for x in v.text.split()])
	forces=np.array(forces)
	return forces
def vasp2xyz():
	try:
		from lxml import etree
	except ImportError:
		print "You need to install python-lxml."
	print "start parse"
	xml = etree.parse("vasprun.xml")
	calculations=xml.xpath('//calculation')
	print 'len(calculations)=',len(calculations);
	from ase.io.trajectory import PickleTrajectory
	atoms=io.read('POSCAR')
	traj = PickleTrajectory('a.traj', 'w',atoms)
	
	print "start find forces and positions"
	allforce=[]
	allpos=[]
	def toarr(v):
		x=v.text.split()
		return map(float,x)
	for i,u in enumerate(calculations):
		print "step : ",i
		forces=map(toarr,u.xpath('./varray[@name="forces"]/v'))
		positions=map(toarr,u.xpath('./structure/varray[@name="positions"]/v'))
		atoms.set_scaled_positions(positions)
		allforce.append(forces)
		allpos.append(positions)
		traj.write()
	np.save('allforce.npy',allforce)
	np.save('allpos.npy',allpos);
	passthru("ase-gui a.traj -o a.xyz ")
def avepos():
	atoms=io.read('POSCAR')
	allpos=np.load('allpos.npy')
	atoms.set_scaled_positions(allpos.mean(axis=0))
	io.write('AVEPOSCAR',atoms,vasp5=True,direct=True,sort=None)
def getmsd():
	allpos=np.load('allpos.npy')
	r=allpos-allpos[0]
	r_atom=np.linalg.norm(r,axis=2)
	r=r_atom.mean(axis=1)
	f=open("msd.txt",'w')
	for i,x in enumerate(r):
		print >>f,"%.3f\t%.3f"%(i*0.25,x)
	with fig('msd.png'):
		pl.plot(np.arange(len(r))*.25,r)
def reducemsd():
	msd=np.loadtxt('600K/msd.txt')
	time=msd[:,0]
	aa=[]
	aa.append(np.loadtxt('msd.txt')[:500,1])
	aa.append(np.loadtxt('600K/msd.txt')[:,1])
	aa.append(np.loadtxt('700K/msd.txt')[:,1])
	aa.append(np.loadtxt('800K/msd.txt')[:,1])
	aa.append(np.loadtxt('900K/msd.txt')[:,1])
	aa=np.array(aa)
	ll=[300,600,700,800,900]
	import matplotlib 
	matplotlib.rcParams['legend.fontsize']=12
	with fig('msd_T.png',legend=True):
		for i,x in enumerate(aa):
			pl.plot(time,x,lw=2,label="%sK"%ll[i])

		pl.xlabel("Time (fs)")
		pl.ylabel("Mean Square Displacement (Angstrom)")
def getcharge():
	atoms=io.read("POSCAR")
	n=len(atoms)+10;
	size=map(int,shell_exec("head -%s CHG|tail -1"%n).split())
	ch=np.loadtxt("CHG",skiprows=n).reshape([size[2],size[1],size[0]])
	ch=np.einsum('kji',ch)
	return ch
	#stress=parseVasprun(vasprun,'stress')
def plotchg():
	ch=getcharge()
	with fig("chg.png"):
		pl.imshow(ch[0])
def getstress():
	
	allpos=np.load('allpos.npy')
	atoms=io.read("POSCAR")
	V=atoms.get_volume()
	n=len(atoms)
	ch=getcharge()
	ch=ch/ch.shape[0]/ch.shape[1]/ch.shape[2]
	#print ch.sum()= total electron number
	print "charge loaded"
	size=ch.shape
	x=np.linspace(0,1,size[0])
	y=np.linspace(0,1,size[1])
	z=np.linspace(0,1,size[2])
	c=itertools.product(x,y,z)
	c=np.array(list(c))
	c=np.einsum('ji',c).reshape([3,size[0],size[1],size[2]])
	r=np.einsum('ijkl,im->mjkl',c,atoms.cell)
	v=np.zeros_like(r)
	stress=[]
	npos=len(allpos)
	from mpi4py import MPI
	comm = MPI.COMM_WORLD
	comm_rank = comm.Get_rank()
	comm_size = comm.Get_size()
	from ase import data
	s=atoms.get_chemical_symbols()
	Z=[data.atomic_numbers[i] for i in s]
	for step in range(comm_rank,npos,comm_size):
		print "rank:",comm_rank
		print "step : ",step
		pos=allpos[step]
		R=np.einsum('ij,jk',pos,atoms.cell)
		for i in range(n):
			print "atoms:",i
			mkdir('./reduce/step%d'%step)
			file='./reduce/step%d/atom%d.npy'%(step,i)
			if exists(file):continue
			for x in range(3):
				v[x]=r[x]-R[i,x]
			p=1/norm(v,axis=0)**3
			q=np.einsum('ijkl,mjkl,jkl->im',v,v,p*ch)
			sigma=np.zeros([3,3])
			for a in range(3):
				for b in range(3):
					s=0.0
					for j in range(n):
						if i==j:continue
						s+=-.5*Z[j]*(R[j,a]-R[i,a])*(R[j,b]-R[i,b])/norm(R[j]-R[i])**3	
					sigma[a,b]=-Z[i]/V*(s+q[a,b])
			np.save(file,sigma)
def gettfc():
	write("""STRUCTURE FILE POSCAR
./POSCAR_unit

FORCE SETS
./FORCE_SETS


SUPERCELL MATRIX PHONOPY
3 0 0
0 3 0
0 0 1""",'input.ph')
	passthru("dynaphopy input.ph OUTCAR --save_force_constants file -r 0.0 7.0 2000 -n 4000")

def getjvq():
	from np.linalg import norm
	unit=io.read("POSCAR_unit")
	atoms=io.read("POSCAR")
	supercell=map(int,norm(atoms.cell,axis=1)/norm(unit.cell,axis=1)+[.5]*3)
	c=supercell
	q=[]
	u=[int(x/2)*2+1 for x in c]
	for i in range(u[0]):
		for j in range(u[1]):
			for k in range(u[2]):
				b=np.array([float(i-c[0]/2)/c[0],float(j-c[1]/2)/c[1],float(k-c[2]/2)/c[2]])
				q.append(b)
	allpos=np.load('allpos.npy')
	v=np.gradient(allpos)[0]
	


def reducestress():
	allpos=np.load('allpos.npy')
	atoms=io.read("POSCAR")
	n=len(atoms)
	m=len(allpos)
	if not exists('sigmas.npy'):
		sigmas=np.zeros([m,n,3,3])
		for step in xrange(m):
			print "step:",step
			for i in range(n):
				file='./reduce/step%d/atom%d.npy'%(step,i)
				sigma=np.load(file)
				sigmas[step,i]=sigma
		np.save('sigmas.npy',sigmas)
	sigmas=np.load('sigmas.npy')
	v=np.gradient(allpos)[0]
	jv=np.einsum('ijkl,ijl->ik',sigmas,v)
	np.save('jv.npy',jv)
	getjvhar()

def getjvhar():
	allpos=np.load('allpos.npy')
	atoms=io.read("POSCAR")
	n=len(atoms)
	m=len(allpos)
	from aces.f import readfc2
	fc2=readfc2()
	sigmas_haR=np.zeros([m,n,3,3])
	R=allpos
	dR=R-np.einsum('i,jk',np.ones(m),atoms.positions)
	dRI=np.einsum('ijk,l',dR,np.ones(n))
	dRJ=np.einsum('ilkj',dRI)
	RI=np.einsum('ijk,l',R,np.ones(n))
	RJ=np.einsum('ilkj',RI)
	sigmas_haR=.5*np.einsum('ijab,tiaj,tibj->tiab',fc2,dRI-dRJ,RI-RJ)
	v=np.gradient(allpos)[0]
	jv_haR=np.einsum('ijkl,ijl->ik',sigmas_haR,v)
	np.save('jv_haR.npy',jv_haR)
def t2c():
	dir1="/home1/xggong/zhouy/tcscripts/bp/nacl.2/0/secondorder/"
	dir2="/home1/xggong/zhouy/tcscripts/bp/nacl.3/0/secondorder/"
	satoms1=io.read(dir1+"SPOSCAR")
	satoms2=io.read(dir2+"SPOSCAR")
	forceset2=np.zeros([2,len(satoms2),3])
	f=open(dir2+"FORCE_SETS")
	for i in range(5):
		f.next()
	for i in range(len(satoms2)):
		force=map(float,f.next().strip().split())
		forceset2[0,i]=force
	for i in range(3):
		f.next()	
	for i in range(len(satoms2)):
		force=map(float,f.next().strip().split())
		forceset2[1,i]=force
	from aces.f import rotationMatrix
	rot=rotationMatrix([1,0,0],-np.pi/4.0)
	rot=rotationMatrix([0,0,1],-np.pi/2.0).dot(rot)
	new_s=satoms2.copy()
	new_s.rotate([1,0,0],-np.pi/4.0,rotate_cell=True)
	new_s.rotate([0,0,1],-np.pi/2.0,rotate_cell=True)
	new_s.write(dir2+'SPOSCAR.1',format='vasp')
def csf():
	"""use compressive sencing to generate force_constant of T
	"""
	from aces.cs1 import runner
	runner(mu=0.0,lam=4.0).run()

def trans_cal():
	from mpi4py import MPI  
	from aces.App import App
	m=App().m
	comm = MPI.COMM_WORLD  
	rank = comm.Get_rank()  
	size = comm.Get_size()  
	print("my rank is: %d" %rank) 
	if rank==0:
		print("Reading force constants from cache")
	d=np.load('fcbin.npz')
	fccenter,fcleft,fcright=d['fccenter'],d['fcleft'],d['fcright']
	
	#fccenter,fcleft,fcright = comm.bcast((fccenter,fcleft,fcright) if rank == 0 else None, root=0)  
	print rank,len(fccenter)
	total=400.0
	fmax=m.fmax
	dm=fmax/total
	intval=dm*size

	omega=np.arange(dm*rank,fmax,intval)#THz
	factor=1e12**2*1e-20*1e-3/1.6e-19/6.23e23
	energies=(omega*2.0*np.pi)**2*factor
	mkdir('tmp')
	from ase.transport.calculators import TransportCalculator
	tcalc =TransportCalculator(h=fccenter,h1=fcleft,h2=fcright,energies=energies,dos=True,logfile='tmp/negf.log'+str(rank))
	if rank==0:
		print ('Calculate Transmission')
	trans=tcalc.get_transmission()
	if rank==0:
		print ('Calculate Dos')
	dos=tcalc.get_dos()*omega
	#np.savez('tmp/result%s.npz'%(rank),x=omega,trans=trans,dos=dos)

	to_txt(['omega','trans','dos'],np.c_[omega,trans,dos],'tmp/result.txt'+str(rank))