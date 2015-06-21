#encoding:utf8
from aces.tools import *
import aces.config as config
from ase.io import read
from ase.io.vasp import write_vasp
from aces.binary import pr
from aces.runners import Runner
class runner(Runner):
	def generate(self):
		m=self.m
		# from minimized structure generate POSCAR
		atoms=read('minimize/range',format='lammps')
		s=atoms.numbers
		symbols=[m.elements[i-1] for i in s ]
		atoms.set_chemical_symbols(symbols)
		write_vasp("POSCAR",atoms,sort="True",direct=True,vasp5=True)
		#generate supercells
		dim=' '.join(str(i) for i in m.supercell)
		passthru(config.phonopy+"-d --dim='%s'"%(dim))
		files=shell_exec("ls *-*").split('\n')
		content="units %s\n"%m.units
		content+="""atom_style      atomic
dimension       3
boundary        p p p 
read_data       structure
%s
%s 
neighbor        2 bin
neigh_modify    every 1 delay 1 check yes
dump 1 all custom 1 dump.force id  fx fy fz xs ys zs
dump_modify 1 format "%%d %%30.20f %%30.20f %%30.20f %%30.20f %%30.20f %%30.20f"
dump_modify  1 sort id
run 0
"""%(m.masses,m.potential)
		from aces.UnitCell.unitcell import UnitCell
		cmd=config.phonopy+"-f "
		maindir=shell_exec('pwd')
		for file in files:
			print file
			dir="dirs/dir_"+file
			mkdir(dir)
			#POSCAR
			mv(file,dir)
			cd(dir)
			#generate structure
			poscar = open(file)
			unit_cell = UnitCell(poscar)
			unit_cell.num_atom_types=len(m.elements)
			write(unit_cell.output_lammps(),"structure")
			#generate in
			write(content,"in")
			#generate dump.force
			shell_exec(config.lammps+" < in >log.out")
			#generate vasprun.xml
			f=open('dump.force')
			for i in range(9):f.next()
			forces=""
			poses=""
			for line in f:
				line=line.split()
				forces+="<v>  %s %s %s </v>\n"%tuple(line[1:4])
				poses+="<v>  %s %s %s  </v>\n"%tuple(line[4:8])
			vasprun='<root><calculation><varray name="forces" >\n'
			vasprun+=forces
			vasprun+='</varray>\n<structure><varray name="positions">\n'+poses
			vasprun+='</varray></structure></calculation></root>\n'
			write(vasprun,'vasprun.xml')
			cmd+=dir+'/vasprun.xml '
			cd(maindir)
		#generate FORCE_SETS
		passthru(cmd)
		#generate mesh.conf
		mesh="""DIM = %s
ATOM_NAME = %s
MP = %s
EIGENVECTORS=.TRUE.
FORCE_CONSTANTS = WRITE
"""%(dim,' '.join(m.elements),' '.join(str(i) for i in m.kpoints))
		write(mesh,'mesh.conf')
		import matplotlib
		matplotlib.use('Agg')
		from matplotlib import pyplot as pl
		passthru(config.phonopy+" --dos  mesh.conf")

		import numpy as np
		xx=np.loadtxt('partial_dos.dat',skiprows=1)
		ndos=len(line)-1
		freq=xx[:,0]
		pdos=xx[:,1:]
		ndos=len(pdos[0,:])
		pl.figure()
		for i in range(ndos):
			pl.plot(freq,pdos[:,i])
		pl.xlabel('Frequency (THz)')
		pl.ylabel('Partial Density of States')
		pl.savefig('partial_dos.png',bbox_inches='tight',transparent=True) 
		pl.figure()
		pl.plot(freq,np.sum(pdos,axis=1),color='red')
		pl.xlabel('Frequency (THz)')
		pl.ylabel('Density of States')
		pl.savefig('total_dos.png',bbox_inches='tight',transparent=True) 
		#calculate paticipation ratio
		
		pr()
		#plot
		xs=[];ys=[]
		for line in open('pr.txt'):
			x,y=map(float,line.split())
			xs.append(x);ys.append(y)
		write("%s"%(sum(ys)/len(ys)),"ave_pr.txt")
		pl.figure()
		pl.plot(xs,ys,'.',color='red')
		pl.ylim([0.0,1.0])
		pl.xlabel('Frequency (THz)')
		pl.ylabel('Paticipation Ratio')
		pl.savefig('Paticipation_atio.png',bbox_inches='tight',transparent=True) 