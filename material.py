# encoding:utf8
# C2N hollow 2D structure
from ase import Atoms,Atom
from math import sqrt,pi
from aces import default
#from ase.io.vasp import write_vasp
from aces.f import writevasp
from ase.io import read
from ase import io
from aces.UnitCell.unitcell import UnitCell
from ase.data import atomic_masses,atomic_numbers
from aces import tools
from aces.modify import get_unique_atoms
from aces.Units import Units
from aces import config
from ase.dft.kpoints import ibz_points
from aces.modify import atoms_from_dump as afd
from aces.lammpsdata import lammpsdata
import numpy as np
from aces.tools import *
from aces.env import SRCHOME,PROJHOME,PROJNAME
class material:
	def __init__(self,opt={}):
		self.__dict__=dict(self.__dict__,**default.default)# all the values needed
		# unit might be changed by opt but need to be used first
		if opt.has_key('units'):
			self.units=opt['units']
		self.units=Units(self.units)
		self.elements=['C','N','B']		
		self.set_parameters()
		self.__dict__=dict(self.__dict__,**opt)
		self.super_setup()
		
	#to be overided	
	def set_parameters(self):
		pass
	def super_setup(self):
		self.units=Units(self.units)
		self.prepare_lammps()

		self.prepare_phonts()
		self.bandpoints=ibz_points['fcc']
		self.bandpoints['X']=[.5,0,0]
		self.bandpoints['Y']=[0,0.5,0]
		self.bandpath=['Gamma','X','Y','Gamma']
		self.premitive=np.eye(3)
		self.dim=self.toString(self.supercell)
		if not self.useS3:
			self.supercell3=self.supercell
		self.setup()
		if self.atomfile:
			atoms=io.read(str(PROJHOME+"/data/"+self.atomfile),format="vasp")
			self.atoms=atoms.repeat([self.latx,self.laty,self.latz])
			self.atoms.center()
		else:
			self.atoms=self.lmp_structure()
		if self.dimension==1:
			self.masses+="\nfix   1d all setforce NULL 0. 0.\nvelocity  all set NULL 0.0 0.0 units box"
		elif self.dimension==2:
			self.masses+="\nfix   1d all setforce NULL NULL 0.\nvelocity  all set NULL NULL 0.0 units box"
	#to be overided
	def setup(self):
		pass
	def prepare_lammps(self):
		self.potential='pair_style	tersoff\npair_coeff	* * %s/BNC.tersoff  %s'%(config.lammpspot,' '.join(self.elements))
		self.dump="dump_modify dump1 element %s"%(' '.join(self.elements))
		masses=self.getMassFromLabel(self.elements)
		self.masses='\n'.join(["mass %d %f"%(i+1,mass) for i,mass in enumerate(masses)])
		m=self
		units=self.units
		m.kb=units.boltz
		m.nktv=units.nktv2p
		if(m.method=="nvt"):m.xp=0;
		m.dtime=m.timestep*100;
		m.tcfactor=units.tcfactor;
		m.excNum=m.aveRate/m.excRate;
		m.swapEnergyRate=m.swapEnergy/(m.excRate*m.timestep);
		
	def prepare_phonts(self):
		masses=self.getMassFromLabel(self.elements)
		self.phontsmasses='\n'.join(["%s %f 0.0"%(label,mass) for label,mass in zip(self.elements,masses)])
		
	def getMassFromLabel(self,labels):
		nums=[atomic_numbers[label] for label in labels]
		masses=[atomic_masses[num] for num in nums]
		return masses
	def toString(self,vec):
		return ' '.join(map(str,vec))
	def extent(self,atoms):
		return atoms.positions.max(axis=0)-atoms.positions.min(axis=0)
		
		
	def structure(self):
		
		self.write()
		
	# to be overrided
	def lmp_structure(self):
		atoms=Atoms()
		return atoms

	# rotate atoms to swap the x z axis for fix=1 and so on, keep the axis right hand
	def swap(self,atoms,fix=1):
		direct=[1,1,1]
		direct[fix]=0
		atoms.rotate(direct,pi,rotate_cell=True)
		order=[[0,2,1],[2,1,0],[1,0,2]][fix]
		cell=atoms.cell[order]
		cell[fix]*=-1
		atoms.set_cell(cell)	
		
	def center_box(self,atoms):
		offset=np.sum(atoms.cell,axis=0)/2
		atoms.translate(-offset)

	def writePOTCAR(self):
		dir='pot'#LDA
		#paw：PAW-LDA
		#paw_gga：PAW-GGA-PW91
		#paw_pbe：PAW-GGA-PBE
		#pot：USPP-LDA
		#pot_GGA：USPP-GGA
		if not self.paw:
			if self.gga:
				dir='pot_GGA'
			else:dir='pot'
		else:
			if not self.gga:
				dir='paw'
			else:
				if self.pbe:
					dir='paw_pbe'
				else:
					dir='paw_gga'
		passthru('cat "" >POTCAR')
		for ele in self.elements:
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
			

	def write(self):
		self.watoms(self.atoms)
	def watoms(self,atoms):
		atoms.write("structure.xyz")
		writevasp(atoms)
		#write_vasp("POSCAR",atoms,sort="True",direct=True,vasp5=True)
		self.POSCAR2data()
		atoms.write('structure.png')	
	def writeatoms(self,atoms,label='atoms'):
		mkcd(label)
		self.watoms(atoms)
		cd('..')
	def getatomicstyle(self):
		a="atom_style atomic"
		if self.creatbonds>0.0:
			a="atom_style bond\natom_modify sort 0 1.\ncomm_modify  cutoff 2.0 "
		return a
	def POSCAR2data(self):
		"""
		poscar = open("POSCAR")
		unit_cell = UnitCell(poscar)
		unit_cell.num_atom_types=len(self.elements)
		tools.write(unit_cell.output_lammps(),"structure")
		"""
		from  ase.io import read
		atoms=read('POSCAR')
		m=self
		atoms.set_pbc([m.xp,m.yp,m.zp])
		#debug(atoms.cell)
		a=lammpsdata(atoms,self.elements)
		rot= a.writedata(filename="structure",creatbonds=self.creatbonds)
		d,p,d1,p1=rot
		np.savetxt('POSCARrot',np.r_[d,p,d1,p1])
		#debug(rot)
		return rot

	def atoms_from_dump(self,filename):
		atoms=afd(filename=filename,elements=self.elements)
		m=self
		atoms.set_pbc([m.xp,m.yp,m.zp])
		return atoms
		
	def dump2POSCAR(self,dumpname,poscar='POSCAR',rotate=True):
		atoms=self.atoms_from_dump(dumpname)
		if rotate:
			rot=np.loadtxt(dirname(dumpname)+'/POSCARrot')
			d,p,d1,p1=rot[:3],rot[3],rot[4:7],rot[7]
			atoms.rotate(d1,-p1,rotate_cell=True)
			atoms.rotate(d,-p,rotate_cell=True)
		#write_vasp(poscar,atoms,sort="True",direct=True,vasp5=True)
		writevasp(atoms,poscar)

	def getboxrange(self):
		file=open("range");
		for i in range(5):
			file.next()
		xlo,xhi=map(float,file.next().split()[:2])
		ylo,yhi=map(float,file.next().split()[:2])
		zlo,zhi=map(float,file.next().split()[:2])
		return (xlo,xhi,ylo,yhi,zlo,zhi);

	def getxrange(self):
		file=open('minimize.xyz');
		n=file.next().split()[0];n=int(n)
		file.next()
		xmin=100000;xmax=-100000;
		ymin=100000;ymax=-100000;
		zmin=100000;zmax=-100000;
		for i in range(n):
			label,x,y,z=file.next().split();
			x,y,z=map(float,[x,y,z]);
			xmin=min(x,xmin);xmax=max(x,xmax);
			ymin=min(y,ymin);ymax=max(y,ymax);
			zmin=min(z,zmin);zmax=max(z,zmax); 
		return (xmin,xmax,ymin,ymax,zmin,zmax);

	def postMini(self):
		xlo,xhi,ylo,yhi,zlo,zhi=self.getboxrange();
		xlo0,xhi0,ylo0,yhi0,zlo0,zhi0=self.getxrange();
		if(self.xp==0):
			xlo=xlo0;xhi=xhi0;
		if(self.yp==0):
			ylo=ylo0;yhi=yhi0;
		if(self.zp==0):
			zlo=zlo0;zhi=zhi0;
		lx=xhi-xlo;ly=yhi-ylo;lz=zhi-zlo;
		if(self.enforceThick):self.zfactor=lz/self.thick;
		else:self.zfactor=1;
		self.S=ly*lz;
		self.box=(xlo,xhi,ylo,yhi,zlo,zhi,lx,ly,lz)

