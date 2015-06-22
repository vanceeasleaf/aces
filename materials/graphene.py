from aces.material import material
from aces.modify import get_unique_atoms
from ase import Atoms,Atom
from math import pi,sqrt
class structure(material):
	def set_parameters(self):
		self.gnrtype='zigzag'
		self.tilt=False
		
	def lmp_structure(self):
		if self.tilt:
			prototype=self.prototype_tilt
		else: prototype=self.prototype
		bond=self.bond
		if self.gnrtype=='armchair':
			col=prototype(self.latx,self.laty)		
		elif self.gnrtype=='zigzag':			
			col=prototype(self.laty,self.latx)	
			col.rotate('z',pi/2,rotate_cell=True)			
		else: raise Exception('Unkown gnr type!')
		col.set_pbc([self.xp,self.yp,self.zp])
		atoms=get_unique_atoms(col)
		cell=atoms.cell*bond
		atoms.set_cell(cell,scale_atoms=True)
		atoms.center()
		return atoms
		
	
	def prototype_tilt(self,latx,laty):
		unit=Atoms()
		unit.append(Atom('C',[1.0/2,0,0]))
		unit.append(Atom('C',[0,sqrt(3)/2,0]))
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
		atoms=Atoms()
		atom=Atoms('C',[(1.0,0.0,0.0)])
		for i in range(6):			
			atom.rotate('z',pi/3*i)
			atoms.extend(atom.copy())
		atoms.set_cell([3.0,sqrt(3),10.0])
		return atoms
			
		