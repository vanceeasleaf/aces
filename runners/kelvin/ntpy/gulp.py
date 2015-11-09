## ## ## gulp.py v.0.1
## ## ## This program executes simple gulp operations.
## ## ## Created: 10/31/2012 - KDP
## ## ## Last Edited: 1/18/2013 - KDP
from os import system
import numpy as np
import strchg as st
import param.const as ct

def freq(strchg, numAtomsUC, gulpName, tempName, gulpExe='gulp'):
	"""
	Executes a gulp program with the string changes neccesary
	and returns the frequencies.

	ntpy.gulp.freq(strchg, numAtomsUC, gulpName, tempName)
	Parameters
	----------
		strchg : dict of type str
			A dictionary holding key:value pairs. Each key
			is original string to be replaced and each value
			is the new string to be subsituted.
		numAtomsUC : int
			The number of atoms in the unit cell.
		gulpName : str
			A string containing the original file name. If
			the file is not included in the pathway then it
			can be the absolute or relative pathway to the
			file.
		tempName : str
			A string containing the new file name. If the
			file is not included in the pathway then it can
			be the absolute or relative pathway to the file.
		gulpExe : str, optional
			A string of the gulp executable. Defaults to 
			"gulp".
	Returns
	----------
		freq : numpy array of type float
			A array returning the frequencies of the k-point.
			The array is of shape (3 * numAtomsUC)
	"""
	# String change the template file
	st.sed(strchg, gulpName, tempName)
	# Execute the new gulp file
	system(gulpExe+ ' <'+ tempName+ ' > output.gulp')
	# system('gulp '+ tempName+ ' output.gulp')
	# Extract the frequencies
	freq = np.zeros( (3 * numAtomsUC), dtype=float)
	freq = np.loadtxt(tempName+ '.freq', comments='--')
	# Remove the .freq file
	system('rm '+ tempName+ '.freq')

	return freq

def eig(strchg, numAtomsUC, gulpName, tempName, kpt, gulpExe='gulp'):
	"""
	Executes a gulp program with the string changes neccesary
	and returns the eigenvectors.

	ntpy.gulp.eig(strchg, numAtomsUC, gulpName, tempName, kpt)
	Parameters
	----------
		strchg : dict of type str
			A dictionary holding key:value pairs. Each key
			is original string to be replaced and each value
			is the new string to be subsituted.
		numAtomsUC : int
			The number of atoms in the unit cell.
		gulpName : str
			A string containing the original file name. If
			the file is not included in the pathway then it
			can be the absolute or relative pathway to the
			file.
		tempName : str
			A string containing the new file name. If the
			file is not included in the pathway then it can
			be the absolute or relative pathway to the file.
		kpt : array of type float
			A numpy array of size 3 that has the three
			dimension k-point to be used.
		gulpExe : str, optional
			A string of the gulp executable. Defaults to 
			"gulp".
	Returns
	----------
		eig : numpy array of type complex
			A eigenvector array where columns correspond to
			modes and rows correspond to atoms in the unit
			cell, where each unit cell atom is split into 
			x, y, z. Shape of array is (3 * numAtomsUC,
			3 * numAtomsUC)
	"""
	numModes = 3 * numAtomsUC
	# String change the template file
	st.sed(strchg, gulpName, tempName)
	# Execute the new gulp file
	system(gulpExe+ ' <'+ tempName+ ' > output.gulp')
	# Grep out eigenvectors
	st.grep(' 1 x', 3 * numAtomsUC, 'output.gulp', 'eigvec_grep.dat')
	xyzDict = dict({ 'x' : ''})
	st.sed(xyzDict, 'eigvec_grep.dat', 'eigvec2.dat')
	system('rm eigvec_grep.dat')
	xyzDict = dict({ 'y' : ''})
	st.sed(xyzDict, 'eigvec2.dat', 'eigvec3.dat')
	system('rm eigvec2.dat')
	xyzDict = dict({ 'z' : ''})
	st.sed(xyzDict, 'eigvec3.dat', 'eigvec4.dat')
	system('rm eigvec3.dat')
	# Read in eigvecs
	if kpt[0] == 0.0 and kpt[1] == 0.0 and kpt[2] == 0.0: 
		# If Gamma Point gulp only prints the real values
		dummy = np.loadtxt('eigvec4.dat', usecols=(1,2,3,4,5,6), comments='--')
		eig = np.zeros( (numModes, numModes), dtype=complex)
		for iSlice in range(0,numModes,3*2):
			eig[:,iSlice:iSlice+(3*2)] = dummy[iSlice*numAtomsUC/2:iSlice*numAtomsUC/2 \
									+(numModes),:]
	else:
		dummy = np.loadtxt('eigvec4.dat', usecols=(1,2,3,4,5,6), comments='--').view(complex)
		eig = np.zeros( (numModes, numModes), dtype=complex)
		for iSlice in range(0,numModes,3):
			eig[:,iSlice:iSlice+3] = dummy[iSlice*numAtomsUC:iSlice*numAtomsUC \
									+numModes,:]
	
	system('rm eigvec4.dat')
	return eig

def vel(strchg, numAtomsUC, gulpName, tempName, kpt, latConst, deltaKpt=1e-3, gulpExe='gulp'):
	"""
	Executes a gulp program with the string changes neccesary
	and returns the velocities.

	ntpy.gulp.vel(strchg, numAtomsUC, gulpName, tempName, deltaKpt=1e-3, gulpExe='gulp')
	Parameters
	----------
		strchg : dict of type str
			A dictionary holding key:value pairs. Each key
			is original string to be replaced and each value
			is the new string to be subsituted.
		numAtomsUC : int
			The number of atoms in the unit cell.
		gulpName : str
			A string containing the original file name. If
			the file is not included in the pathway then it
			can be the absolute or relative pathway to the
			file.
		tempName : str
			A string containing the new file name. If the
			file is not included in the pathway then it can
			be the absolute or relative pathway to the file.
		kpt : array of type float
			A numpy array of size 3 that has the three
			dimension k-point to be used.
		latConst : float
			A float of the lattice constant in angstroms.
		deltaKpt : float, optional
			A float that determines what difference to use
			for the difference theorem. Defaults to 10e-5.
		gulpExe : str, optional
			A string of the gulp executable. Defaults to 
			"gulp".
	Returns
	----------
		vel : numpy array of type float
			Returns the velocity for the kpt in a array of
			shape (3 * numAtomsUC, 3)
	"""
	def _kptstrchg(strchg, kpt):
		strchg['KPT'] = '{:.10f} {:.10f} {:.10f}'.format(kpt[0], kpt[1], kpt[2])
		return strchg

	vel = np.zeros( (3 * numAtomsUC, 3), dtype=float)

	# Convert to m/s
	velConv = ((1.0 / ct.value('centi')) * ct.value('c')) / (1.0 / (latConst * ct.value('ang')))

	# For all three directions
	for idim in range(3):
		if kpt[idim] == 0.5: # kpt at right boundary
			freqVal = freq(strchg, numAtomsUC, gulpName, tempName, gulpExe=gulpExe) 
			# Change kpt
			kpt[idim] = kpt[idim] - deltaKpt
			strchg = _kptstrchg(strchg, kpt)
			freqValMdk = freq(strchg, numAtomsUC, gulpName, tempName, gulpExe=gulpExe)
			vel[:, idim] = ((freqVal - freqValMdk) / deltaKpt) * velConv
			# Reset kpt
			kpt[idim] = kpt[idim] + deltaKpt
			strchg = _kptstrchg(strchg, kpt)
		elif kpt[idim] == -0.5: # kpt at left boundary
			freqVal = freq(strchg, numAtomsUC, gulpName, tempName, gulpExe=gulpExe)
			# Change kpt
			kpt[idim] = kpt[idim] + deltaKpt
			strchg = _kptstrchg(strchg, kpt)
			freqValPdk = freq(strchg, numAtomsUC, gulpName, tempName, gulpExe=gulpExe)
			vel[:, idim] = ((freqValPdk - freqVal) / deltaKpt) * velConv
			# Reset kpt
			kpt[idim] = kpt[idim] - deltaKpt
			strchg = _kptstrchg(strchg, kpt)
		elif kpt[idim] == 0.0: # kpt at gamma point
			freqVal = freq(strchg, numAtomsUC, gulpName, tempName, gulpExe=gulpExe)
			# Change kpt
			kpt[idim] = kpt[idim] + deltaKpt
			strchg = _kptstrchg(strchg, kpt)
			freqValPdk = freq(strchg, numAtomsUC, gulpName, tempName, gulpExe=gulpExe) 
			vel[:, idim] = ((freqValPdk - freqVal) / deltaKpt) * velConv
			# Reset kpt
			kpt[idim] = kpt[idim] - deltaKpt
			strchg = _kptstrchg(strchg, kpt)
		else:
			freqVal = freq(strchg, numAtomsUC, gulpName, tempName, gulpExe=gulpExe) 
			# Change kpt
			kpt[idim] = kpt[idim] + deltaKpt
			strchg = _kptstrchg(strchg, kpt)
			freqValPdk = freq(strchg, numAtomsUC, gulpName, tempName, gulpExe=gulpExe) 
			kpt[idim] = kpt[idim] - (2.0 * deltaKpt)
			strchg = _kptstrchg(strchg, kpt)
			freqValMdk = freq(strchg, numAtomsUC, gulpName, tempName, gulpExe=gulpExe) 
			vel[:, idim] = ((freqValPdk - freqValMdk) / (2.0 * deltaKpt)) * velConv
			# Reset kpt
			kpt[idim] = kpt[idim] + deltaKpt
			strchg = _kptstrchg(strchg, kpt)

	return vel


