## ## ## Kevin Parrish - gulp run script ## ## ##
import nmd
import numpy as np
import os.path
from os import system
import sys
sys.path.append(nmd.ntpyPath) # Needed to recognize ntpy module
import ntpy.strchg as st
import ntpy.gulp as gp

def removeFile(filename):
	if os.path.isfile(filename):
		system('rm '+ filename)

def nmdGulp(latVec, latType, dim, mass, freqConv, gulpName, gulpTrans, \
			gulpExe, numModes, numAtomsUC, kpt):

	freq = np.zeros( (kpt[:,0].size, numModes), dtype=float)
	eigvec = np.zeros( (kpt[:,0].size, numModes, numModes), dtype=complex)
	vel = np.zeros( (kpt[:,0].size, numModes, 3), dtype=float)

	for ikpt in range(kpt[:,0].size):
		strchg = dict({ 'tempName' : gulpTrans,
					'MASS' : '{0:10.10f}'.format(mass),
					'KPT' : '{0:.10f} {1:.10f} {2:.10f}'.format(kpt[ikpt,0], \
							kpt[ikpt,1], kpt[ikpt,2]),
					'ALAT' : '{0:.10f} {1:.10f} {2:.10f}'.format(latVec, \
							latVec, latVec)
					})
	
		freq[ikpt,:] = gp.freq(strchg, numAtomsUC, gulpName, gulpTrans, \
								gulpExe=gulpExe) * freqConv

		vel[ikpt,:,:] = gp.vel(strchg, numAtomsUC, gulpName, gulpTrans, \
								kpt[ikpt,:], latVec, gulpExe=gulpExe)

		eigvec[ikpt,:,:] = gp.eig(strchg, numAtomsUC, gulpName, gulpTrans, \
								kpt[ikpt,:], gulpExe=gulpExe)

	return freq, vel, eigvec
#-- END nmdGulp --#


###--- MAIN ---###

## ## GULP ## ##
skip = False
if len(sys.argv) > 1:
	assert sys.argv[1] == '--nogulp', 'Invalid option \"'+ str(sys.argv[1])+ '\"'
	print "Gulp run skipped"
	skip = True
if not skip:	
	print "Begin gulp run"
	freq, vel, eigvec = nmdGulp(nmd.realLat, nmd.latType, nmd.dim, nmd.realMass, nmd.freqConv, \
						nmd.gulpName, nmd.gulpTrans, nmd.gulpExe, nmd.numModes, nmd.numAtomsUC, \
						nmd.kpt)

	eigvec = np.reshape(eigvec, (nmd.numKpts * nmd.numModes, nmd.numModes))

	filename = 'post.gulp.npz'
	np.savez(filename, freq=freq, vel=vel, eigvec=eigvec)
# Cleanup
removeFile('input.gulp')
removeFile('output.gulp')
removeFile(nmd.gulpTrans+ '.freq')


## ## LAMMPS ## ##
print "Begin lammps file creation"

# File for single submission script
removeFile(nmd.lammpsSubRunName+ '.sh')
lmpF = open(nmd.lammpsSubRunName+ '.sh', 'a')

for iseed in range(nmd.numSeeds):
	ext = '.'+ str(iseed)

	## Extra Lammps String Subsitutions
	nmd.lmpInFile['MAIN_LOG_FILE'] = 'out.lammps.main'+ ext
	nmd.lmpInFile['IN_POS'] = nmd.inPosName
	nmd.lmpInFile['LMP_TMP'] = nmd.lammpsRunName+ext
	nmd.lmpInFile['SEED'] = str((iseed + 1) * 11111)
	nmd.lmpInFile['FFT_LOG_FILE'] = 'out.lammps.fft'+ ext
	nmd.lmpInFile['OUT_VEL'] = 'out.lammps.vel'+ ext
	nmd.lmpInFile['W_STEP'] = nmd.wStep
	nmd.lmpInFile['T_FFT'] = nmd.tFft
	nmd.lmpInFile['T_TOTAL'] = nmd.tTotal

	lmpSubFile = dict({ 'NUM_CPU' : nmd.lmpCpu,
					'LMP_TEMP' : nmd.lammpsRunName+ ext,
					'EXECPATH' : nmd.lammpsExecPath,
					'EXECUTABLE' : 'lmp_generic',
					'RUNPATH' : nmd.runpath
					})

	# Lammps Run file
	st.sed(nmd.lmpInFile, nmd.lammpsName, nmd.lammpsRunName+ ext)

	# Lammps Submission file
	st.sed(lmpSubFile, nmd.lammpsSubName, nmd.lammpsSubRunName+ ext+ '.sh')

	# Make a single submission script
	lmpF.write(r'qsub -l walltime='+ str(nmd.lmpWallTime)+ ':00:00 -l nodes=1:ppn='+ \
				str(nmd.lmpCpu)+ ' '+ str(nmd.lammpsSubRunName)+ ext+ '.sh\n')
lmpF.close()


## ## NMD CALCULATION ## ##
print "Begin nmd file creation"

# File for single submission script
removeFile(nmd.nmdSubRunName+ '.sh')
nmdF = open(nmd.nmdSubRunName+ '.sh', 'a')

for iseed in range(nmd.numSeeds):
	for ikslice in range(nmd.numKslice):
		seed = '.'+ str(iseed)
		ext = seed+ '.'+ str(ikslice)
		## NMD Processing Sub files
		nmdSubFile = dict({ 'NUM_CPU' : nmd.pyCpu,
							'LMP_TEMP' : nmd.nmdRunName+ ext+ '.py',
							'EXECPATH' : nmd.pythonExecPath,
							'EXECUTABLE' : 'python',
							'RUNPATH' : nmd.runpath
							})
		st.sed(nmdSubFile, nmd.nmdSubName, nmd.nmdSubRunName+ ext+ '.sh')

		## NMD Processing Run files
		nmdRunFile = dict({ 'SEED' : seed,
							'IKSLICE' : ikslice,
							})
		st.sed(nmdRunFile, nmd.nmdName, nmd.nmdRunName+ ext+ '.py')

		# Make a single submission script
		nmdF.write(r'qsub -l walltime='+ str(nmd.pyWallTime)+ ':00:00 -l nodes=1:ppn='+ \
					str(nmd.pyCpu)+ ',mem='+ str(nmd.pyMem)+ 'gb '+ str(nmd.nmdSubRunName)+ ext+ '.sh\n')
nmdF.close()

