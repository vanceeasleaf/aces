from ase import io
import numpy as np 
atoms=io.read('minimize/POSCAR')
bulk=np.abs(atoms.positions[:,0]-20)<20 
knot=np.abs(atoms.positions[:,0]-70)<20
index=np.arange(len(atoms),dtype='int')-10
def filt(a):
	return a[(a>=0) *(a<520)]
from aces.dos import plot_regiondos
plot_regiondos([(filt(index[bulk]),'bulk'),(filt(index[knot]),'knot')])