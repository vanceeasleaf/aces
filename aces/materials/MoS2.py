from aces.materials  import Material
from aces.modify import get_unique_atoms
from ase import Atoms,Atom
from math import pi,sqrt
from ase.dft.kpoints import ibz_points
class structure(Material):
	def set_parameters(self):
		self.gnrtype='zigzag'
		self.tilt=False
		self.bond=3.193/sqrt(3)
	def setup(self):
		self.elements=["Mo","S"]
		self.bandpoints=ibz_points['hexagonal']
		self.bandpath=['Gamma','M','K','Gamma']
	def lmp_structure(self):
		if self.tilt:
			prototype=self.prototype_tilt
		else: prototype=self.prototype
		bond=self.bond
		if self.gnrtype=='armchair':
			col=prototype(self.latx,self.laty)		
		elif self.gnrtype=='zigzag':			
			col=prototype(self.laty,self.latx)		
			self.swap(col,2)
		else: raise Exception('Unkown gnr type!')	
		col.center()	
		atoms=get_unique_atoms(col)
		col.set_pbc([self.xp,self.yp,self.zp])
		cell=atoms.cell*bond
		atoms.set_cell(cell,scale_atoms=True)
		atoms.center()

		return atoms
		
	
	def prototype_tilt(self,latx,laty):
		unit=Atoms()
		b=3.134/self.bond/2
		unit.append(Atom('Mo',[1.0/2,0,0]))
		unit.append(Atom('S',[0,sqrt(3)/2,b]))
		unit.append(Atom('S',[0,sqrt(3)/2,-b]))
		unit.set_cell((3.0/2,sqrt(3),10.0))
		unit.cell[0,1]=sqrt(3)/2
		col=unit.repeat((latx,laty,1))
		return col
	
	def prototype(self,latx,laty):
		#armchair
		unit=self.ring()
		col=unit.repeat((latx,laty,1))
		return col
		
	def ring(self):
		#armchair ring
		b=3.134/self.bond/2
		atoms=Atoms()
		atom=Atoms('MoS2',[(1.0,0.0,0.0),(-1.0,0,b),(-1.0,0,-b)])
		for i in range(3):			
			atom.rotate('z',pi*2.0/3*i)
			atoms.extend(atom.copy())
		atoms.set_cell([3.0,sqrt(3),10.0])
		return atoms
			
		