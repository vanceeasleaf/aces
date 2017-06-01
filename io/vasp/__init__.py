# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2017-06-01 21:49:49
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-01 21:52:16
import numpy as np
from aces.f import toString
def writevasp(atoms,file='POSCAR'):
	f=open(file,'w')
	s=np.array(atoms.get_chemical_symbols())
	ss=atoms.get_scaled_positions()
	print >>f,'ACES POSCAR'
	print >>f,'1.0'
	for x in atoms.cell:
		print >>f,toString(x)
	#ele=np.unique(s)
	ele=[]
	for a in s:
		if a in ele:
			continue
		ele.append(a)
	print >>f,toString(ele)
	a=[]
	p=np.arange(len(s))
	for e in ele:
		a.append(p[s==e])
	ns=[len(x) for x in a]
	print >>f,toString(ns)
	print >>f,'Direct'
	v=[]
	for x in a:
		for u in x:
			v.append(u)
			print >>f,toString(ss[u])
	f.close()
	x=np.array(v,dtype=np.int).argsort()
	np.savetxt('POSCARswap',x)