# -*- coding: utf-8 -*-
# @Author: YangZhou
# @Date:   2015-10-14 20:43:32
# @Last Modified by:   YangZhou
# @Last Modified time: 2016-12-02 12:56:05
from aces.tools import *
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
	from ase import io
	from ase.io.trajectory import PickleTrajectory
	atoms=io.read('POSCAR')
	traj = PickleTrajectory('a.traj', 'w',atoms)
	
	print "start find forces and positions"
	def toarr(v):
		x=v.text.split()
		return map(float,x)
	for i,u in enumerate(calculations):
		print "step : ",i
		forces=map(toarr,u.xpath('./varray[@name="forces"]/v'))
		positions=map(toarr,u.xpath('./structure/varray[@name="positions"]/v'))
		atoms.set_scaled_positions(positions)
		traj.write()
	
	passthru("ase-gui a.traj -o a.xyz ")

	#stress=parseVasprun(vasprun,'stress')