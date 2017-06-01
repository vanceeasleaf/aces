import numpy as np
from aces.tools import parseyaml
def pr():
    data =parseyaml('mesh.yaml')
    file=open("pr1.txt",'w')
    for phonon in data['phonon']:
        for band in phonon['band']:
            frequency=band['frequency']
            eigenvector=np.array(band['eigenvector'])
            u=ny.vectorize(lambda x:np.sum(x*x),eigenvector)
            pr=1.0/(u*u).sum()
            print >>file,"%f\t%f"%(frequency,pr)
    file.close()
if __name__=='__main__':
    pr();


