#encoding:utf8
from aces.tools import *
import aces.config as config
from ase.io import read
from ase.io.vasp import write_vasp
from aces.binary import pr
from aces.runners import Runner
from aces.graph import plot,series
import numpy as np
from aces.runners.phonopy import Runner as Prunner
class runner(Runner):
	def generate(self):
		m=self.m
		Prunner(m).run()
		
		

