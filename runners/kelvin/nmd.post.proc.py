## ## ## Kevin Parrish - 7/26/2012	## ## ##
## ## ## Post Processing for nmd	## ## ##

import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('/home/kevin/projects/ntpy')
# from ntpy.leastsqbound import leastsqbound as lsb
import ntpy.param.const as ct
# import ntpy.param.lj as lj

# Add path
repoPath = '/home/kevin/projects/isostrain/nmd/'
sys.path.append(repoPath) # Needed to import nmd

### Functions ###

def bteRTAtau(heatCap, vel, life, V):
	return vel * vel * life * heatCap / V

def meanFreePath(life, vel):
	magVel = np.sqrt(vel[:,:,0]**2 + vel[:,:,1]**2 + vel[:,:,2]**2)
	return life * magVel

def accumFunc(life, vel, heatCap, V):
	# Find mean free path
	mfp = meanFreePath(life, vel)
	# Find conductivities
	condx = bteRTAtau(heatCap, vel[:,:,0], life, V)
	
	# Reshape to 1d arrays
	mfp1d = np.reshape(mfp, -1)
	condx1d = np.reshape(condx, -1)

	# Sort by mean free path
	mfpOrd = np.lexsort((condx1d[:], mfp1d[:]))

	return np.cumsum(condx1d[mfpOrd]), mfp1d[mfpOrd]

### NMD Processing ### 
print "Begin NMD processing"


job = ['20K.0.4x', '20K.0.6x', '20K.0.8x']
rSizes = np.array( [4, 6, 8] )
rstrains = [0.0]

## Extrapolate Thermal Conductivities ##
tmpKappa = np.zeros( (1, 1, 3) )
kappa = np.zeros( (1, 1) )

for itemp in range(1):
	for istrain in range(1):
		for isize in range(3):
			# Import run data
			file = np.load('post.nmd.'+ job[isize]+ '.npz')
			tmpKappa[itemp, istrain, isize] = file['condx']

		values = np.polyfit(1.0 / rSizes, 1.0 / tmpKappa[itemp,istrain,:], 1)
		kappa[itemp, istrain] = 1.0 / values[1]
		print kappa[itemp, istrain]
		np.savez('post.nmd.cond.npz', kappa=kappa, tmpKappa=tmpKappa, strains=rstrains)


## Plot accumulation function
for itemp in range(1):
	for istrain in range(1):
		for isize in range(3):
			# Import run data
			file = np.load('post.nmd.'+ job[isize]+ '.npz')
			life = file['life']
			vel = file['vel']
			V = file['V']
			heatCap = ct.value('kb')
			print 'life\t'+ str(np.max(life[:,:]))
			print 'vel\t'+ str(np.max(vel[:,:,:]))
			print 'heatCap\t'+ str(heatCap)
			print 'V\t'+ str(V)
			accumCond, accumMfp = accumFunc(life, vel, heatCap, V)

			plt.semilogx(accumMfp, accumCond) 
			plt.xlim(xmin=1e-9)
			plt.title('NMD Accumulation Function')
			plt.xlabel(r'$\Lambda \/ (m)$')
			plt.ylabel(r'$k \/ \left(\frac{W}{m-K}\right)$')
			plt.savefig('pic.accum.'+ job[isize]+ '.png')
			plt.show()
			plt.clf()

