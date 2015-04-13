# encoding:utf8
# C2N hollow 2D structure
from ase import Atoms,Atom
from math import sqrt,pi
from aces import default
from ase.io.vasp import write_vasp
from aces.UnitCell.unitcell import UnitCell

class structure:
	def __init__(self,home,opt):
		self.home=home
		self.__dict__=dict(self.__dict__,**default.default)# all the values needed
		self.potential='pair_style        tersoff\npair_coeff      * * %s/potentials/BNC.tersoff  C N'%home
		self.masses="mass 1 12.01\nmass 2 14.00"
		self.dump="dump_modify dump1 element C N"
		self.latysmall=3; # 窄边包括2*latysmall+1排C原子，latysmall=3时，为7 AGNR*/
		self.latExtend=5;# 相对于窄边，长边的y方向两头都延伸latExtend个六边形，latExtend=2时，长边为[(latysamll-1)+latExtend]*2+1=13 AGNR*/
		self.latxsmall=7;#奇数
		self.latxbig=12;#最简单模型，中间是长边，两头是等长的短边, latxsmall表示短边由几条*折线*组成，偶数*/
		self.yp=0
		self.fixud=1
		self.__dict__=dict(self.__dict__,**opt)
	def extent(self,atoms):
		    xmin=100000;xmax=-100000;
		    ymin=100000;ymax=-100000;
		    zmin=100000;zmax=-100000;
		    for pos in atoms.positions:
		        x,y,z=pos
		        xmin=min(x,xmin);xmax=max(x,xmax);
		        ymin=min(y,ymin);ymax=max(y,ymax);
		        zmin=min(z,zmin);zmax=max(z,zmax); 
		        lx=xmax-xmin;
		        ly=ymax-ymin;
		    return (lx,ly);
	def structure(self):
		latysmall,latExtend,latxsmall,latxbig,bond=[int(self.latysmall),int(self.latExtend),int(self.latxsmall),int(self.latxbig),float(self.bond)]
		if(latxsmall%2==0):latxsmall+=1;
		if(latxbig%2==1):latxbig+=1;
		atoms=self.agnr(latysmall,latxsmall,0);
		unit=self.agnr(latysmall-1+2*latExtend,latxbig,0)
		unit.translate([latxsmall*1.5,-(latExtend-0.5)*sqrt(3),0])
		atoms.extend(unit)
		unit=self.agnr(latysmall,latxsmall,1);
		unit.translate([(latxsmall+latxbig)*1.5,0,0])
		atoms.extend(unit)
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
		unit_cell.num_atom_types=2
		lammps=open("structure","w")
		lammps.write(unit_cell.output_lammps())
		lammps.close()
		#data("structure").write("new")
