## ## ## Kevin Parrish ## ## ##
import nmd
import numpy as np

###--- MAIN ---###
loadSpecED = np.zeros( (nmd.numSeeds, nmd.numKslice, nmd.sliceLength, nmd.numTstep/2, \
					nmd.numModes), dtype=float )

preSpecED = np.zeros( (nmd.numSeeds, nmd.numKpts, nmd.numTstep/2, nmd.numModes), dtype=float)

avgSpecED = np.zeros( (nmd.numKpts, nmd.numTstep/2, nmd.numModes), dtype=float)

specED = np.zeros( (nmd.irrKpt[:,0].size, nmd.numTstep/2, nmd.numModes), dtype=float )

## Load the specED
for iseed in range(nmd.numSeeds):
	for ikslice in range(nmd.numKslice):
		ext = '.'+ str(iseed)+ '.'+ str(ikslice)
		# Load nmd files for averaging
		nmdFile = np.load('out.nmd'+ ext+ '.npy')
		loadSpecED[iseed,ikslice,:,:,:] = nmdFile
	#-- END ikslice --#
#-- END iseed --#

## Consolidate the k-point slices
for ikslice in range(nmd.numKslice):
	currSlice = ikslice * nmd.sliceLength # Used for indexing to the current slice
	preSpecED[:,currSlice:currSlice+nmd.sliceLength,:,:] = loadSpecED[:,ikslice,:,:,:]
#-- END ikslice --#

## Average over seeds
for iseed in range(nmd.numSeeds):
	avgSpecED[:,:,:] = avgSpecED[:,:,:] + preSpecED[iseed,:,:,:]
#-- END iseed --#
avgSpecED[:,:,:] = avgSpecED[:,:,:] / nmd.numSeeds

## Use cubic symmetries to average
for ikpt in range(nmd.numKpts):
	specED[nmd.degen[ikpt],:,:] = specED[nmd.degen[ikpt],:,:] + avgSpecED[ikpt,:,:]
#-- END ikpt --#

for irr in range(nmd.irrKpt[:,0].size):
	specED[irr,:,:] = specED[irr,:,:] / nmd.irrKptCnt[irr]
#-- END irr --#

np.save('post.nmd.npy', specED)

