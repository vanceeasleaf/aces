#encoding:utf8
import sys
from aces.Units import Units
from aces.tools import *
import aces.config as config
from ase.io import read
from ase.io.vasp import write_vasp
from aces.runners import Runner
import numpy as np
class runner(Runner):
	def runcmd(self):
		return config.mpirun+" %s "%self.m.cores+config.phonts+"  >log.out"
	def generate(self):
		m=self.m
		coordination=self.phontsAtoms()
		content0="species %d\n"%(len(m.elements))+m.phontsmasses+"""
D3_cutoff 2.0
kpoints %s 1
delta 0.005
numerical_2der T
numerical_3der T
iter_steps 3
pdos 0. 70. 200 10.
temp 60. 400. 10
AbInitio  T F 
FP_interface LAMMPS
#phonons_only T
Lattice  1.0
%s
end
"""%(m.toString(m.kpoints),coordination)
		write(content0,'phonons_input.dat')
		passthru(config.phonts) # generate many displacement files
		mkdir('lammps');cd('lammps')
		content="units %s\n"%m.units
		content+="""atom_style      charge
dimension       3
boundary        p p p 
read_data       GENERIC
%s
%s 
neighbor        1.1 bin
neigh_modify    every 1 delay 1 check yes
dump 1 all custom 1 *.dump id  fx fy fz 
dump_modify 1 format "%%d %%30.20f %%30.20f %%30.20f"
dump_modify  1 sort id
run 0
"""%(m.masses,m.potential)
		shell_exec("mv ../*.str .")
		strs=shell_exec("ls *.str").split("\n")
		for str in strs:
			dir=str.replace("str","dir")
			mkdir(dir)
			write(content.replace("GENERIC",str),dir+"/in")
			mv(str,"%s/%s"%(dir,str))
			cd(dir)
			passthru(config.lammps+" <in >out.dat")
			cd('..')
		dirs=shell_exec("ls |grep dir").split("\n")
		for dir in dirs:
			print dir
			cp(dir+"/0.dump","../"+dir.replace("dir","out"))
		cp("1.0000.dir/out.dat","../1.0000.out")
		cd('..')
		content0=content0.replace('AbInitio  T F','AbInitio  F T')
		write(content0,'phonons_input.dat')
		passthru(config.phonts)

	def phontsAtoms(self):
		m=self.m
		m.dump2POSCAR('minimize/range')
		atoms=read('POSCAR')
		cell=atoms.get_cell()
		if not np.allclose(np.diag(np.diag(cell)),cell):
			raise Exception('phonts needs cell to be orthorgnal')
		content="cell %f %f %f\n"%tuple(np.diag(cell))
		content+="natoms %d\n"%(len(atoms))
		content+="fractional\n"
		pos=atoms.get_scaled_positions()
		for i,atom in enumerate(atoms):
			content+="%s %s\n"%(atom.symbol,m.toString(pos[i]))
		return content		