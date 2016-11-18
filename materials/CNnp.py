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
class structure:
	#姨妈巾结构
	def __init__(self,opt):
		self.__dict__=dict(self.__dict__,**default.default)# all the values needed
		self.set_parameters()
		self.__dict__=dict(self.__dict__,**opt)
		self.setup()
		
	def set_parameters(self):
		self.elements=['C','N']
		self.latysmall=3; # 窄边包括2*latysmall+1排C原子，latysmall=3时，为7 AGNR*/
		self.latExtend=5;# 相对于窄边，长边的y方向两头都延伸latExtend个六边形，latExtend=2时，长边为[(latysamll-1)+2*latExtend]*2+1=13 AGNR*/
		self.latxsmall=7;#奇数
		self.latxbig=12;#最简单模型，中间是长边，两头是等长的短边, latxsmall表示短边由几条*折线*组成，偶数*/
		self.yp=0
		self.fixud=0	
		self.seg=0
	def setup(self):
		self.potential='pair_style        tersoff\npair_coeff      * * %s/potentials/BNC.tersoff  %s'%(SRCHOME,' '.join(self.elements))
		self.masses=""
		self.phontsmasses=""
		i=1
		for a in self.elements:
			num=atomic_numbers[a]
			mass=atomic_masses[num]
			self.masses+="mass %d %f\n"%(i,mass)
			self.phontsmasses+="%s %f 0.0\n"%(a,mass)
			i+=1
		self.dump="dump_modify dump1 element %s"%(' '.join(self.elements))
		
	def extent(self,atoms):
		xmax=atoms.positions[:,0].max()
		xmin=atoms.positions[:,0].min()
		ymax=atoms.positions[:,1].max()
		ymin=atoms.positions[:,1].min()
		lx=xmax-xmin;
		ly=ymax-ymin;
		return (lx,ly);
			
	def structure(self):
		latysmall,latExtend,latxsmall,latxbig,bond=[int(self.latysmall),int(self.latExtend),int(self.latxsmall),int(self.latxbig),float(self.bond)]
		if(latxsmall%2==0):latxsmall+=1;
		if(latxbig%2==1):latxbig+=1;
		atoms=self.agnr(latysmall,latxsmall+1,0);
		unit=self.agnr(latysmall-1+2*latExtend,latxbig,0)
		unit.translate([latxsmall*1.5,-(latExtend-0.5)*sqrt(3),0])
		atoms.extend(unit)
		unit=self.agnr(latysmall,latxsmall+1,0);
		unit.translate([(latxsmall+latxbig-1)*1.5,0,0])
		atoms.extend(unit)
		temp=Atoms()	
		for i in range(self.seg):
			unit=atoms.copy()
			unit.translate([(latxsmall+latxbig-1)*1.5*i,0,0])
			temp.extend(unit)
		atoms.extend(temp)
		atoms=get_unique_atoms(atoms)
		lx,ly=self.extent(atoms)
		atoms.set_cell([lx,ly,100])
		atoms.set_cell([lx*bond,ly*bond,100],scale_atoms=True)
		atoms.set_pbc([1,1,1])
		atoms.center(vacuum=10*bond)
		self.atoms=atoms
		
		self.write()
		print 'read_data structure'
	

	"""
	* 生成宽度为n类型为type的折线的坐标,并把它放在第p列
	* @author zhouy
	* @input n,type,p
	* @output 包含2*n+1个原子的数组
	"""
	def zhexian(self,n,p,type):
		m=n*2+1;
		atoms=Atoms()
		for i in range(m):
			y=i*sqrt(3)/2
			x=0.5*(1+(type*2-1)*(i%2)-type)+p*1.5
			atoms.append(Atom('C',(x,y,0.0)))
		
		return atoms
	
	"""
	* 生成宽度为n，长度为p的AGNR的坐标,第一列折线的类型由type决定
	* @author zhouy
	* @input n,p,type
	* @output (2*n+1)*p个原子的数组
	"""
	def agnr(self,n,p,type):
		atoms=Atoms()
		for i in range(p):
			unit=self.zhexian(n,i,type+(1-2*type)*(i%2));
			atoms.extend(unit)			
		return atoms

	def write(self):
		self.atoms.write("CN.xyz")
		write_vasp("POSCAR",self.atoms,sort="True",direct=True,vasp5=True)
		poscar = open("POSCAR")
		unit_cell = UnitCell(poscar)
		unit_cell.num_atom_types=len(self.elements)
		tools.write(unit_cell.output_lammps(),"structure")
		#data("structure").write("new")

