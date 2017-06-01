# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2017-06-01 21:53:01
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-01 21:57:40
from aces.f import *
def writefc2(fc2,filename='FORCE_CONSTANTS_2ND'):
	natom=len(fc2)
	s="%d\n"%natom
	for i in range(natom):
		for j in range(natom):
			s+="%d\t%d\n"%(i+1,j+1)
			s+=matrixFormat(fc2[i,j])
	write(s,filename)
def writefc3(fc3,atoms,satoms,filename='FORCE_CONSTANTS_3RD'):
	psm=premitiveSuperMapper(atoms,satoms)
	natom=len(fc3)
	s=""
	n=0
	scale=np.amax(np.abs(fc3))
	for i in range(natom):
		for j in range(natom):
			for k in range(natom):
				if np.allclose(np.zeros([3,3,3]),fc3[i,j,k]/scale,atol=1e-04):continue
				n+=1
				s+="\n%d\n"%n
				s+=psm.maps2p(i,j,k)
				s+=matrix3Format(fc3[i,j,k])
	s="%d\n"%n+s			
	write(s,filename)
class premitiveSuperMapper:
	def __init__(self,atoms,satoms):
		self.atoms=atoms
		self.satoms=satoms
	def getS2p(self):
		atoms=self.atoms
		cellp=atoms.cell
		posp=atoms.get_scaled_positions()
		satoms=self.satoms
		cell=satoms.cell
		pos=satoms.get_scaled_positions()

		#which cell is i j in
		vpos=pos*np.linalg.norm(cell,axis=1)/np.linalg.norm(cellp,axis=1)
		v=np.floor(vpos)
		
		celloffset=v.dot(cellp)
		vpos-=v
		s2p=-np.ones(len(vpos),dtype='int')
		for pd,p in enumerate(vpos):
			for id,a in enumerate(posp):
				if np.allclose(p,a):
					s2p[pd]=id
					break

		return s2p,v,celloffset
	def maps2p(self,i,j,k):
		if not hasattr(self, 's2p'):
			self.s2p,c,self.celloffset=self.getS2p()
		v=self.celloffset

		s=toString(v[j]-v[i])+'\n'
		s+=toString(v[k]-v[i])+'\n'
		s+='%d %d %d\n'%(self.s2p[i]+1,self.s2p[j]+1,self.s2p[k]+1)
		return s
def matrixFormat(mat):
	n,m=mat.shape
	s=""
	for k in range(n):
		for p in range(m):
			s+="\t%f"%mat[k,p]
		s+="\n"	
	return s
def matrix3Format(mat):
	n,m,q=mat.shape
	s=""
	for k in range(n):
		for p in range(m):
			for t in range(q):
				s+=" %d  %d  %d       %f\n"%(k+1,p+1,t+1,mat[k,p,t])
	return s
def refinefc3():
	f=open('FORCE_CONSTANTS_3RD')
	g=open('fc3new','w')
	nblock=int(f.next().split()[0])
	fc3=[]
	r1t=[]
	r2t=[]
	idx=[]
	#print >>g,nblock
	for i in range(nblock):
		f.next()
		f.next()
		r1=np.array(map(float,f.next().split()))
		r1t.append(r1)
		r2=np.array(map(float,f.next().split()))
		r2t.append(r2)
		idx.append(map(int,f.next().split()))
		fc=np.zeros([3,3,3])
		for i in range(3):
			for j in range(3):
				for k in range(3):
					fc[i,j,k]=float(f.next().split()[3])
		fc3.append(fc)
	scale=np.amax(np.abs(fc3))
	u=[]
	for i in range(nblock):
		if np.allclose(np.zeros([3,3,3]),fc3[i]/scale,atol=1e-01*.3) :
			u.append(False)
		else:
			u.append(True)
	u=np.array(u)
	fc3=np.array(fc3)[u]
	r1t=np.array(r1t)[u]
	r2t=np.array(r2t)[u]
	idx=np.array(idx)[u]
	v=np.unique(idx[:,0])
	filters=np.array([idx[:,0]==c for c in v])
	
	for i,f in enumerate(filters):
		a=np.arange(len(idx))[f][0]
		fc3[a]-= fc3[f].sum(axis=0)	
	v=np.unique(idx[:,1])
	filters=np.array([idx[:,1]==c for c in v])
	
	for i,f in enumerate(filters):
		a=np.arange(len(idx))[f][0]
		fc3[a]-= fc3[f].sum(axis=0)	
	nblock=len(fc3)
	print >>g,nblock
	for i in range(nblock):
		print >>g,i+1
		print >>g,toString(np.around(r1t[i],3))
		print >>g,toString(np.around(r2t[i],3))
		print >>g,toString(idx[i])
		print >>g,matrix3Format(fc3[i])
def rotatefc2(t,direct=[0,0,1],file1='FORCE_CONSTANTS_2ND',file2='fc2new'):
	fc2=readfc2(file1)
	fc=rotate_fc2(fc2,t,direct) 
	writefc2(fc,file2)
def rotate_fc2(fc2,t,direct=[0,0,1]):
	M=rotationMatrix(direct,t*np.pi/180.0)
	fc=np.einsum('im,klmn->klin',M,fc2)
	fc=np.einsum('in,klmn->klmi',M,fc)
	return fc
def rotatefc3(t,direct=[0,0,1],file1='FORCE_CONSTANTS_3RD',file2='fc3new'):
	f=open(file1)
	g=open(file2,'w')
	M=rotationMatrix(direct,t*np.pi/180.0)
	nblock=int(f.next().split()[0])
	print >>g,nblock
	print >>g,""
	for i in range(nblock):
		f.next()
		print >>g,f.next(),
		r1=np.array(map(float,f.next().split()))
		r1t=M.dot(r1)
		print >>g,toString(r1t)
		r1=np.array(map(float,f.next().split()))
		r1t=M.dot(r1)
		print >>g,toString(r1t)
		print >>g,f.next(),
		fc=np.zeros([3,3,3])
		for i in range(3):
			for j in range(3):
				for k in range(3):
					fc[i,j,k]=float(f.next().split()[3])

		fc=np.einsum('ij,jkl->ikl',M,fc)
		fc=np.einsum('ij,kjl->kil',M,fc)
		fc=np.einsum('ij,klj->kli',M,fc)
		print >>g,matrix3Format(fc)
def read_forces(filename):
	try:
	    from lxml import etree as ElementTree
	    xmllib="lxml.etree"
	except ImportError:
	    try:
	        import xml.etree.cElementTree as ElementTree
	        xmllib="cElementTree"
	    except ImportError:
	        import xml.etree.ElementTree as ElementTree
	        xmllib="ElementTree"

    """
    Read a set of forces on atoms from filename, presumably in
    vasprun.xml format.
    """
    calculation=ElementTree.parse(filename
                                  ).getroot().find("calculation")
    for a in calculation.findall("varray"):
        if a.attrib["name"]=="forces":
            break
    nruter=[]
    for i in a.getchildren():
        nruter.append([float(j) for j in i.text.split()])
    nruter=np.array(nruter,dtype=np.double)
    return nruter
def readfc2(filename='FORCE_CONSTANTS'):
	f=open(filename)
	line=f.next()
	natom=int(line)
	fc=np.zeros([natom,natom,3,3])
	for i in range(natom):
		for j in range(natom):
			f.next()
			for k in range(3):
				fc[i,j,k]=map(float,f.next().split())
	return fc

def readfc3(atoms,unit,filename='FORCE_CONSTANTS_3RD'):
	f=open(filename)
	nblock=int(f.next().split()[0])
	n=len(atoms)
	fc3=np.zeros([n,n,n,3,3,3])
	for i in range(nblock):
		f.next() #blank
		f.next(), #index
		r1=np.array(map(float,f.next().split()))
		x=r1/norm(unit.cell,axis=1)
		r2=np.array(map(float,f.next().split()))
		y=r2/norm(unit.cell,axis=1)
		u=.5*(x-np.abs(x))+.5*(y-np.abs(y))
		idx=np.array(map(int,f.next().split()))-np.array([1,1,1])
		#print idx
		pos1=-u.dot(unit.cell)+unit.positions[idx[0]]
		pos2=(x-u).dot(unit.cell)+unit.positions[idx[1]]
		pos3=(y-u).dot(unit.cell)+unit.positions[idx[2]]
		ii=hash(pos1,atoms.positions)
		jj=hash(pos2,atoms.positions)
		kk=hash(pos3,atoms.positions)
		fc=np.zeros([3,3,3])
		for i in range(3):
			for j in range(3):
				for k in range(3):
					fc[i,j,k]=float(f.next().split()[3])
		fc3[ii,jj,kk]=fc
	return fc3
def hash(pos,positions):
	for i,pos in enumerate(positions):
		if(np.allclose(pos,positions[i])):
			return i
