#encoding: utf-8 
#a class to calculate inequality of atoms in a given structure
import sys
import math
from ase import Atoms
from ase.io import read
from ase.calculators.neighborlist import NeighborList
from numpy import dot
from math import atan2,pi,floor
#sys.path.append('/home1/xggong/zhouy/tcscript/pizza');
#from dump import dump

		
class inequality:
	def __init__(self):
		pass
	def readDump(self):
		atoms=read('../range',format='lammps')
		atoms. set_pbc((1, 1, 1))
		self.atoms=atoms	
		self.natom=len(atoms)
			
	def buildNeighborList(self):
		atoms=self.atoms
		nl = NeighborList([0.8 for atom in atoms],self_interaction=False,bothways=True)
		nl.update(atoms)
		self.nl=nl
		neigh=[]
		for i in xrange(self.natom):
			neigh.append(self.sort_nei(i))
			#print neigh[i]
		self.neigh=neigh
			
	def nonequ5(self):		

		# 不等价的cluster个数
		nonequ=[];
		atoms=self.atoms
		for i in xrange(16):
			nonequ.append(0)
		for i in xrange(self.natom):
			code=0;
			for j in xrange(len(self.neigh[i])):
				k=self.neigh[i][j];
				dex=1
				if atoms.numbers[k]==1:dex=0
				code+=(1<<j)*dex
			dex=1
			if atoms.numbers[i]==1:dex=0
			code+=(1<<3)*dex
			nonequ[code]+=1;
		return max(nonequ)/(self.natom+0.0);	
		
	def run(self):
		self.readDump()
		self.buildNeighborList()
		print self.nonequ5()	
		
		
	def sort_nei(self,i):
		atoms=self.atoms
		nl=self.nl
		occu=[0,0,0,0];
		indices, offsets = nl. get_neighbors(i)
		for j, offset in zip(indices, offsets):
			pos= atoms. positions[j] + dot(offset, atoms. get_cell())-atoms.positions[i]
			ang=atan2(pos[1],pos[0])+pi
			idx=int(ang/(pi/2.0));
			
			# 原子j占据了第N象限那个位置*/
			occu[idx]=j;
			
		if(occu[2]==0):
			return [occu[0],occu[1],occu[3]];
		else:
			return [occu[2],occu[3],occu[1]]	

				
if __name__=='__main__':
	ie=inequality()
	ie.run()
