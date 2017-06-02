# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2017-06-01 21:49:49
# @Last Modified by:   YangZhou
# @Last Modified time: 2017-06-02 21:31:17
import numpy as np
from aces.f import toString
from aces.tools import *
from aces import config
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

def writePOTCAR(options,elements):
	dir='pot'#LDA
	#paw：PAW-LDA
	#paw_gga：PAW-GGA-PW91
	#paw_pbe：PAW-GGA-PBE
	#pot：USPP-LDA
	#pot_GGA：USPP-GGA
	if not options.paw:
		if options.gga:
			dir='pot_GGA'
		else:dir='pot'
	else:
		if not options.gga:
			dir='paw'
		else:
			if options.pbe:
				dir='paw_pbe'
			else:
				dir='paw_gga'
	passthru('cat "" >POTCAR')
	for ele in elements:
		file=config.vasppot+"/%s/%s/POTCAR"%(dir,ele)
		z=False
		if not exists(file):
			file+='.Z'
			z=True
		assert exists(file)
		if z:
			passthru('zcat %s >> POTCAR'%file)
		else:
			passthru('cat %s >> POTCAR'%file)
	#s=''.join([tools.read(config.vasppot+"/%s/%s/POTCAR.Z"%(dir,ele)) for ele in self.elements])
	#tools.write(s,'POTCAR')
def parseVasprun(vasprun,tag="forces"):
	collection = []
	for event, element in vasprun:
			if element.attrib['name'] == tag:
	 			for v in element.xpath('./v'):
	 				collection.append([float(x) for x in v.text.split()])
	collection=np.array(collection)
	return collection
def writeKPOINTS(kpoints):
		s="""A
0
Monkhorst-Pack
%s
0  0  0
	"""%' '.join(map(str,kpoints))
		write(s,'KPOINTS')