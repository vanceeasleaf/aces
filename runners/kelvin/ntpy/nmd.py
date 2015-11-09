## ## ## nmd.py v.0.1
## ## ## This program performs nmd specific functions
## ## ## Created: 02/21/2013 - KDP

import struct
import numpy as np
import numpy.ma as ma

def lmpReadBin(filename, numTstep, numAtoms, numCol=3, oldStyle=False, bteOrd='@'):
	"""
	This function reads in a lammps binary file.

	ntpy.lmpReadBin(filename, numTstep, numAtoms, numCol=3, oldStyle=False, bteOrd='@')
	Parameters
	----------
		filename : str
			The name of the lammps binary file.
		numTstep : int
			How many timesteps to read in.
		numAtoms : int
			How many atoms in the dump file.
		numCol : int, optional
			The number of parameters that were dumped.
		oldStyle : bool, optional
			If True, then the binary is assumed to follow the old style. That is,
			the ntimestep and natoms variables are not of ctype int64_t and the
			boundary variables are omitted. Otherwise, the current binary format
			is assumed (default).
		bteOrd : {'@','<','>'}, optional
			Designates the byte order (endianness). Character '<' reads using 
			little-endian while '>' reads using big-endian. By default, '@' is
			used which follows the machine endianness.
	Returns
	-------
		values : list of ndarrays of type float
			A list of shape (numCol) with numpy arrays of shape (numTstep, numAtoms).
	"""
	# Initialize qdot
	if oldStyle:
		bIt = 'q'
		bIb = 8
		boundary = True
	else:
		bIt = 'i'
		bIb = 4
		boundary = False
		
	values = [0] * numCol
	for dims in range(numCol):
		values[dims] = np.zeros( (numTstep, numAtoms), dtype=float)

	f = open(filename, 'rb')

	# Loop over snapshots in file
	for itime in range(numTstep):
		real = f.read(bIb)
		if real:
			ntimestep = struct.unpack(bteOrd+ bIt, real)[0]
			natoms = struct.unpack(bteOrd+ bIt, f.read(bIb))[0]
			triclinic = struct.unpack(bteOrd+ 'i', f.read(4))[0]
			if boundary is True:
				boundary0 = struct.unpack(bteOrd+ 'i', f.read(4))[0]
				boundary1 = struct.unpack(bteOrd+ 'i', f.read(4))[0]
				boundary2 = struct.unpack(bteOrd+ 'i', f.read(4))[0]
				boundary3 = struct.unpack(bteOrd+ 'i', f.read(4))[0]
				boundary4 = struct.unpack(bteOrd+ 'i', f.read(4))[0]
				boundary5 = struct.unpack(bteOrd+ 'i', f.read(4))[0]
			xlo = struct.unpack(bteOrd+ 'd', f.read(8))[0]
			xhi = struct.unpack(bteOrd+ 'd', f.read(8))[0]
			ylo = struct.unpack(bteOrd+ 'd', f.read(8))[0]
			yhi = struct.unpack(bteOrd+ 'd', f.read(8))[0]
			zlo = struct.unpack(bteOrd+ 'd', f.read(8))[0]
			zhi = struct.unpack(bteOrd+ 'd', f.read(8))[0]
			if triclinic:
				xy = struct.unpack(bteOrd+ 'd', f.read(8))[0]
				xz = struct.unpack(bteOrd+ 'd', f.read(8))[0]
				yz = struct.unpack(bteOrd+ 'd', f.read(8))[0]
			size_one = struct.unpack(bteOrd+ 'i', f.read(4))[0]
			nchunk = struct.unpack(bteOrd+ 'i', f.read(4))[0]

			# Loop over processor chunks in file
			for i in range(nchunk):
				n = struct.unpack(bteOrd+ 'i', f.read(4))[0]
			
			# Check if number of values matches with user arguments
			assert n == (numAtoms * numCol), 'ERROR: numAtoms or numCol does not agree with file'

			# Read chunk
			for i in range(n / numCol):
				for dims in range(numCol):
					values[dims][itime,i] = struct.unpack(bteOrd+ 'd', f.read(8))[0]
		else:
			print 'ERROR: EOF reached'
			return 1
	###--- END itime ---###

	return values
###--- END lmpReadBin ---###


def nmdProc(velx, vely, velz, eigVec, latPosx, latPosy, latPosz,
			latVecx, latVecy, latVecz, kpt, ikpt, imode, numAtomsUC, numUC, mass,
			numTstep):
	"""
	This function performs normal mode decomposition.

	nmd.nmdProc(velx, vely, velz, eigVec, latPosx, latPosy, latPosz, latVecx, latVecy,
				latVecz, kpt, ikpt, imode, numAtomsUC, numUC, mass, numTstep)
	Parameters
	----------
		velx : ndarray of type float
			Numpy array containing the x velocities.
		vely : ndarray of type float
			Numpy array containing the y velocities.
		velz : ndarray of type float
			Numpy array containing the z velocities.
		eigVec : ndarray of type float
			Numpy array containing the eigenvectors. Eigenvectors must be in shape
			(numKpts * numModes, numModes), where the numModes by numModes arrays
			are stacked ontop of each other.
		latPosx : ndarray of type float
			Numpy array containing the x positions.
		latPosy : ndarray of type float
			Numpy array containing the y positions.
		latPosz : ndarray of type float
			Numpy array containing the z positions.
		latVecx : float
			The lattice vector in the x direction.
		latVecy : float
			The lattice vector in the y direction.
		latVecy : float
			The lattice vector in the z direction.
		kpt : ndarray of type float
			The array of kpoints of shape (numKpts, 3).
		ikpt : int
			The index of the current kpoint in kpt.
		imode : int
			The index of the current mode.
		numAtomsUC : int
			The number of atoms in the unit cell.
		numUC : int
			The number of unit cells.
		mass : float or ndarray
			Either the mass of all the atoms or a numpy array of the masses of the
			unit cell or the entire simulation.
		numTstep : int
			The number of times steps.
	Returns
	-------
		specEDFft : ndarray
			Numpy array containind the spectral energy density of the fft and of 
			shape (numTstep / 2) and type float.
	"""
	# Initialize qdot
	qdot = np.zeros( (numTstep) )

	# Prereference functions
	tile = np.tile
	conjugate = ma.conjugate

	# Spatial fourier transform factor
	spatial = 2.0 * np.pi * 1j * (\
		latPosx[:]*( (kpt[ikpt,0])/(latVecx) ) + \
		latPosy[:]*( (kpt[ikpt,1])/(latVecy) ) + \
		latPosz[:]*( (kpt[ikpt,2])/(latVecz) ) ) 
	
	# Conjugate eigenvectors
	eigx = tile(conjugate(eigVec[(numAtomsUC*3*ikpt)+0: \
				(numAtomsUC*3*(ikpt+1)):3, imode]),numUC)
	eigy = tile(conjugate(eigVec[(numAtomsUC*3*ikpt)+1: \
				(numAtomsUC*3*(ikpt+1)):3, imode]),numUC)
	eigz = tile(conjugate(eigVec[(numAtomsUC*3*ikpt)+2: \
				(numAtomsUC*3*(ikpt+1)):3, imode]),numUC)

	# qdot: normal mode kinetic energy corrdinate
	qdot = np.sum(((velx * eigx) + (vely * eigy) + (velz * eigz)) * \
			np.exp(spatial) * np.sqrt(mass/numUC), axis=1)

	# keXcorr: kinetic energy autocorrelation
	result = np.correlate(qdot, qdot, mode='full', old_behavior=False)
	keXcorr = result[result.size/2:] / result[result.size/2]

	# keFft: kinetic energy FFT
	keFft = np.fft.fft(keXcorr[:])

	# specEDFft: spectral energy density for a single FFT
	specEDFft = (keFft[:].real * keFft[:].real) + \
			(keFft[:].imag * keFft[:].imag)

	# Add spectral energy denisty for a single FFT to whole
	return specEDFft[:numTstep/2]

#-- END nmdProc --#
