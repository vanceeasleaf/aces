#encoding:utf8
from aces.runners import Runner
from ase.io.vasp import write_vasp
from ase import io
from aces.tools import *
import numpy as np
from  aces import config
from aces.runners.phonopy import runner as Runner
from aces.runners.cs import read_forces
def matrixFormat(mat):
	n,m=mat.shape
	s=""
	for k in range(n):
		for p in range(m):
			s+="\t%f"%mat[k,p]
		s+="\n"	
	return s
class runner(Runner):
	def generate(self):
		#self.minimizePOSCAR()
		shell_exec('cp minimize/POSCAR .')
		self.get_almin()
		self.displacements()
		files,self.NDATA=self.getfiles()
		self.getvasprun(files)
		self.getdispforce(files)
		self.get_fitin()
		self.getfcs()
		self.get_anphonin()
		self.run_anphonin()
	def run_anphonin(self):
		passthru(config.anphon+"< band.in > band.out")
		passthru(config.anphon+"< dos.in > dos.out")
		passthru(config.anphon+"< tc.in > tc.out")
	def getfcs(self):
		passthru(config.alm+"< fit.in > fit.out")
		assert exists("alm.fcs")
		assert exists("alm.xml")
	def getfiles(self):
		files=shell_exec("ls *.POSCAR|sort -n").split('\n')
		assert len(files)>0 and not files[0]==''
		return files,len(files)
	def get_fitin(self):
		content=read('alm.in').replace('suggest','fitting')
		fitting="""&fitting
	NDATA = %d
	DFILE = disp_all.dat
	FFILE = force_all.dat
/
"""%self.NDATA
		write(content+fitting,'fit.in')
	def getdispforce(self,files):
		force=""
		disp=""
		orig=io.read('POSCAR-supercell').positions
		for dir0 in files:
			forcearr=read_forces('dirs/dir_%s/vasprun.xml'%dir0)/25.7110# in Rd/bohr
			force+=matrixFormat(forcearr)
			disparr=(io.read('dirs/dir_%s/POSCAR'%dir0).positions-orig)*1.889726# in bohr
			disp+=matrixFormat(disparr)
		write(disp,"disp_all.dat")
		write(force,"force_all.dat")
	def displacements(self):
		passthru(config.alm+"< alm.in > alm.out")
		files=shell_exec("ls *pattern*").split()
		passthru(config.almdisp+self.m.toString(files))

	def get_almin(self):
		m=self.m
		m.atoms=io.read('POSCAR')
		atoms=m.atoms.repeat(m.supercell)
		general="""&general
  PREFIX = alm
  MODE = suggest
  NAT = %s; NKD = %s
  KD = %s
/
"""%(len(atoms),len(m.elements),m.toString(m.elements))
		interaction="""&interaction
  NORDER = 1
/
"""
		cell="""&cell
  1.889726
  %s
/
"""%('\n  '.join([m.toString(atoms.cell[i]) for i in range(3)]))
		cutoff="""&cutoff
  *-* None 
/
"""
		pos='  \n'.join(['%s '%(m.elements.index(a.symbol)+1)+m.toString(atoms.get_scaled_positions()[i]) for i,a in enumerate(atoms)])
		position="""&position
	%s
/
"""%pos
		write(general+interaction+cell+cutoff+position,'alm.in')
		write_vasp('POSCAR-supercell',atoms,sort="True",direct=True,vasp5=True)
	def get_anphonin(self):
		m=self.m
		masses=m.getMassFromLabel(m.elements)
		general="""&general
  PREFIX = alm
  MODE = phonons
  FCSXML = alm.xml
  NKD = %s
  KD = %s
  MASS = %s
/
"""%(len(m.elements),m.toString(m.elements),m.toString(masses))
		cell="""&cell
	1.889726
  %s
/
"""%('\n\t'.join([m.toString(m.atoms.cell[i]) for i in range(3)]))
		bp=m.bandpoints
		s=""
		for i in range(len(m.bandpath)-1):
			x1,x2=m.bandpath[i],m.bandpath[i+1]
			s+="  %s "%x1[0]+m.toString(bp[x1])+" %s "%x2[0]+m.toString(bp[x2])+" 101\n"
		kpoint="""&kpoint
  1
%s/
"""%s
		write(general+cell+kpoint,'band.in')
		kpoint="""&kpoint
  2
  %s
/
"""%m.toString(m.kpoints)
		write(general+cell+kpoint,'dos.in')
		write(general.replace('phonons','RTA')+cell+kpoint,'tc.in')