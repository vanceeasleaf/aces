#encoding:utf8
from aces.runners import Runner
from aces.runners.correlation import runner as Crun
from aces.runners.phonopy import runner as Prun
import numpy as np
import phonopy.file_IO as file_IO
import dynaphopy.functions.iofunctions as reading
import dynaphopy.classes.controller as controller
class runner(Runner):
	def generate(self):
		crun=Crun(self.m)
		print "running phonopy to generate FORCE_SETS"
		prun=Prun(self.m)
		prun.run()
		structure = reading.read_from_file_structure_poscar('POSCAR')
		structure.set_force_set(file_IO.parse_FORCE_SETS(filename='FORCE_SETS'))
		structure.set_primitive_matrix([[0.5, 0.0, 0.0],
                                [0.0, 0.5, 0.0],
                                [0.0, 0.0, 0.5]])
		structure.set_super_cell_phonon(np.diag(self.m.supercell))
		trajectory = reading.read_from_file_trajectory('/home/abel/VASP/Si-dynamic_300/RUN1/OUTCAR',structure)
		calculation = controller.Calculation(trajectory)
		
		crun.run()
		crun.vd.life_yaml(correlation_supercell=self.m.correlation_supercell)