#encoding:utf8
from aces.tools import *
import aces.config as config
from ase.io import read
from ase.io.vasp import write_vasp
from aces.binary import pr
from aces.runners import Runner
from aces.UnitCell.unitcell import UnitCell
from aces.graph import plot,series
from aces.script.vasprun import exe as lammpsvasprun
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pl
import time
import numpy as np
from aces.runners.phonopy import runner as Runner
class runner(Runner):
	def generate(self):
		m=self.m
		cp('minimize/POSCAR','.')
		self.getVaspRun_vasp()
	def q(self):
		a=shell_exec("grep TOTEN OUTCAR |tail -1").split("=")[1].strip().replace("eV","")
		print self.m.ecut,a