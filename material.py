# encoding:utf8
# C2N hollow 2D structure
from ase import Atoms,Atom
from math import sqrt,pi
from aces import default
from ase.io.vasp import write_vasp

from aces.UnitCell.unitcell import UnitCell
from ase.data import atomic_masses,atomic_numbers
from aces import tools
from aces.modify import get_unique_atoms
from aces.env import SRCHOME
from aces.Units import Units
from aces import config
class material:
	def __init__(self,opt):
		self.__dict__=dict(self.__dict__,**default.default)# all the values needed
		self.elements=['C','N','B']
		self.set_parameters()
		self.__dict__=dict(self.__dict__,**opt)
		self.super_setup()
		
	#to be overided	
	def set_parameters(self):
		pass
	def super_setup(self):
		self.prepare_lammps()
		self.prepare_phonts()
		self.setup()
	#to be overided
	def setup(self):
		pass
	def prepare_lammps(self):
		self.potential='pair_style	tersoff\npair_coeff	* * %s/potentials/BNC.tersoff  %s'%(SRCHOME,' '.join(self.elements))
		self.dump="dump_modify dump1 element %s"%(' '.join(self.elements))
		masses=self.getMassFromLabel(self.elements)
		self.masses='\n'.join(["mass %d %f"%(i+1,mass) for i,mass in enumerate(masses)])
		m=self
		units=Units(m.units)
		m.kb=units.boltz
		m.nktv=units.nktv2p
		if(m.method=="nvt"):m.xp=0;
		m.dtime=m.timestep*100;
		m.tcfactor=units.tcfactor;
		m.excNum=m.aveRate/m.excRate;
		m.swapEnergyRate=m.swapEnergy/(m.excRate*m.timestep);
		m.units=units	
	def prepare_phonts(self):
		masses=self.getMassFromLabel(self.elements)
		self.phontsmasses='\n'.join(["%s %f 0.0"%(label,mass) for label,mass in zip(self.elements,masses)])
		
	def getMassFromLabel(self,labels):
		nums=[atomic_numbers[label] for label in labels]
		masses=[atomic_masses[num] for num in nums]
		return masses
		
	def extent(self,atoms):
		xmax=atoms.positions[:,0].max()
		xmin=atoms.positions[:,0].min()
		ymax=atoms.positions[:,1].max()
		ymin=atoms.positions[:,1].min()
		lx=xmax-xmin;
		ly=ymax-ymin;
		return (lx,ly);
		
	def structure(self):
		self.atoms=self.lmp_structure()
		self.write()
		
	# to be overrided
	def lmp_structure(self):
		atoms=Atoms()
		return atoms
		
		
		
	def writePOTCAR(self):
		s=''.join([tools.read(config.vasppot+"/%s/POTCAR"%ele) for ele in self.elements])
		tools.write(s,'POTCAR')
			

	def write(self):
		self.atoms.write("CN.xyz")
		write_vasp("POSCAR",self.atoms,sort="True",direct=True,vasp5=True)
		self.POSCAR2data()
	
	def POSCAR2data(self):
		"""
		poscar = open("POSCAR")
		unit_cell = UnitCell(poscar)
		unit_cell.num_atom_types=len(self.elements)
		tools.write(unit_cell.output_lammps(),"structure")
		"""
		from aces.lammpsdata import lammpsdata
		from  ase.io import read
		atoms=read('POSCAR')
		a=lammpsdata(atoms,self.elements)
		a.writedata()
		
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
