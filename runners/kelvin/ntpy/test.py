import struct
import numpy as np
import nmd

def lmpReadBin(filename, numTstep, numAtoms, numDim=3):
	values = [0] * numDim
	for dims in range(numDim):
		values[dims] = np.zeros( (numTstep, numAtoms), dtype=float)

	f = open(filename, 'rb')

	# Loop over snapshots in file
	for itime in range(numTstep):
		real = f.read(8)
		if real:
			ntimestep = struct.unpack('q', real)[0]
			natoms = struct.unpack('q', f.read(8))[0]
			triclinic = struct.unpack('i', f.read(4))[0]
			boundary0 = struct.unpack('i', f.read(4))[0]
			boundary1 = struct.unpack('i', f.read(4))[0]
			boundary2 = struct.unpack('i', f.read(4))[0]
			boundary3 = struct.unpack('i', f.read(4))[0]
			boundary4 = struct.unpack('i', f.read(4))[0]
			boundary5 = struct.unpack('i', f.read(4))[0]
			xlo = struct.unpack('d', f.read(8))[0]
			xhi = struct.unpack('d', f.read(8))[0]
			ylo = struct.unpack('d', f.read(8))[0]
			yhi = struct.unpack('d', f.read(8))[0]
			zlo = struct.unpack('d', f.read(8))[0]
			zhi = struct.unpack('d', f.read(8))[0]
			if triclinic:
				xy = struct.unpack('d', f.read(8))[0]
				xz = struct.unpack('d', f.read(8))[0]
				yz = struct.unpack('d', f.read(8))[0]
			size_one = struct.unpack('i', f.read(4))[0]
			nchunk = struct.unpack('i', f.read(4))[0]

			# Loop over processor chunks in file
			for i in range(nchunk):
				n = struct.unpack('i', f.read(4))[0]
			
			# Check if number of values matches with user arguments
			assert n == (numAtoms * numDim), 'ERROR: numAtoms or numDim does not agree with file'

			# Read chunk
			for i in range(n / numDim):
				for dims in range(numDim):
					values[dims][itime,i] = struct.unpack('d', f.read(8))[0]
		else:
			print 'ERROR: EOF reached'
			return 1
	###--- END itime ---###
	return values
###--- END lmpReadBin ---###
					
			

a, b, c = readVel('out.lammps.vel.0.1.bin', nmd.numTstep, nmd.numAtoms)
# print a
print a.shape
print b.shape
print c.shape
print a
print b
print c
