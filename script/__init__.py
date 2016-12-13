# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2015-10-14 20:43:32
# @Last Modified by:   YangZhou
# @Last Modified time: 2016-12-13 14:54:35
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
	jv_haR=np.einsum('ijkl,ijl->ik',sigmas_haR,v)
	np.save('jv_haR.npy',jv_haR)

