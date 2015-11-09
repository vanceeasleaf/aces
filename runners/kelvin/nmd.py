## ## ## Kevin Parrish		## ## ##
## ## ## nmd.py config file	## ## ##
import numpy as np
import sys
# ntpyPath = '/opt/mcgaugheygroup/ntpy'
ntpyPath = '/home/kevin/projects/ntpy'
sys.path.append(ntpyPath) # Needed to recognize ntpy module
import ntpy.lattice as lt
import ntpy.param.lj as lj
import ntpy.param.const as ct

## ## ## BEGIN USER INPUT	## ## ##

## Lattice Parameters
dim = [4, 4, 4]
latType = 'fcc'
realMass = 1.0 * lj.value('mass') * ct.value('kilo') * ct.value('avog')
ljMass = 1.0
realLat = lj.value('lat20') * (1.0+ 0.0)
ljLat = ((lj.value('lat20') / lj.value('sigma')) * 1e-10) * (1.0 + 0.0)

## Name Parameters

gulpName = 'gulp.template'				# Gulp Template Name
gulpTrans = 'input.gulp'				# Gulp Transient Name
gulpExe = 'gulp_ifort'					# Gulp Executable Name

lammpsName = 'lammps.template'			# Lammps Template Name
lammpsRunName = 'in.lammps'				# Lammps Run Name

lammpsSubName = 'lammps.sub.template'	# Lammps Submission Template Name
lammpsSubRunName = 'tmp.lammps.sub'		# Lammps Submission Run Name

lammpsExecPath = '/opt/mcgaugheygroup/lammps-2Nov10/src'	# Lammps Executable Path

grepName = 'lammps.vel.grep'			# Grep Template Name
grepRunName = 'in.grep'					# Grep Run Name

grepSubName = lammpsSubName				# Grep Submission Template Name
grepSubRunName = 'tmp.grep.sub'			# Grep Submission Run Name

nmdName = 'nmd.proc.template'			# NMD Template Name
nmdRunName = 'in.nmd.proc'				# NMD Run Name

nmdSubName = lammpsSubName				# NMD Submission Template Name
nmdSubRunName = 'tmp.nmd.sub'			# NMD Submission Run Name

pythonExecPath = '/opt/python/python-2.7.3/bin'	# Python Executable Path

inPosName = 'in.pos'					# Lammps Position File Name

## Time Specifications
tTotal = 2**21							# Total Run Time
tFft = 2**19							# Time for each FFT Block (Sets longest freq)
wStep = 2**5							# Time Step for Sampling (Sets shortest freq)
dt = 0.002								# Time Step for Lammps

## Run Specifications (used for queue submission)
pyWallTime = 2							# NMD Python Processing Run Time
pyCpu = 1								# NMD Python Number of CPUs
pyMem = 2								# NMD Python Memory Requirement
lmpWallTime = 2							# Lammps Run Time
lmpCpu = 1								# Lammps Number of CPUs
numSeeds = 10							# Number of seeds to use
numKslice = 4							# Number of k-point Slices (used to reduce memory)
freqConv = 2.0 * np.pi * ct.value('c') \
	* ct.value('tocenti') * lj.value('tau')	# Frequency conversion for gulp

## Lammps String Subsitutions (more included in create.files.py)
lmpInFile = dict({ 'W_STEP' : wStep,
				'T_FFT' : tFft,
				'T_TOTAL' : tTotal,
				'TEMP_PARAM' : str(20),
				'LAT_PARAM' : str(ljLat * dim[0]),
				'STRAIN_PARAM' : str(0.0)
				})


## ## ## END USER INPUT		## ## ##
## ## ## NB: USER SHOULD	## ## ##
## ## ## NOT EDIT BELOW!!!	## ## ##

## Name Parameters (cont.)
runpath = str(sys.path[0])				# Current path

## Time Specifications (cont.)
numFft = tTotal / tFft					# Number of FFT Blocks
numTstep = tFft / wStep					# Number of FFT Time Steps

## Freqency Specifications (used for plotting)
freqStep = 2 * np.pi / (tFft * dt)		# Frequency Step
freqMax = 2 * np.pi / (wStep * dt * 2)	# Maximum Frequency
numFreq = tFft / (wStep * 2)			# Number of Frequencies

omega = np.zeros( (numFreq) )
omega = np.arange(numFreq) * (freqMax / numFreq)

## Create lammps position file ##
mylamp = lt.Lammps(inPosName)
myBox = lt.Block([ljLat, ljLat, ljLat], latType, dim, atom_mass = [ljMass])
mylamp.buildLammps([myBox])
latPos = lt.buildNumpy([myBox])

## Lattice properties
numAtomsUC = myBox.numAtomsUcell()
numUcell = myBox.numUcell()
numAtoms = myBox.numAtoms()
Lx = myBox.Lx()
Ly = myBox.Ly()
Lz = myBox.Lz()
numModes = 3 * numAtomsUC

## Create kpoints
kpt = lt.kpt(dim, conv=True)
kptIndex = np.zeros( (kpt[:,0].size) )
kptIndex = range(kpt[:,0].size)

# Other kpoint parameters
numKpts = kpt[:,0].size
sliceLength = kpt[:,0].size / numKslice

# Create the irreducible k-point list
irrKpt, irrKptCnt, degen = lt.cubicKptSym(kpt)

