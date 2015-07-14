root="/home/xggong/home1/zhouy/"
def exepath(a,abs=False):
	if abs:
		return ' '+a+' '
	else:
		return ' '+dirpath(a)+' '
def dirpath(a):
	return root+a.lstrip('/')

php       = exepath("php/bin/php")
lammps    = exepath("lammps-24Apr13/src/lmp_ubuntu")
lammpspot = dirpath("lammps-24Apr13/potentials")
mpirun    = exepath("/opt/intel/mpi/openmpi/1.6.3/icc.ifort/bin/mpirun  -np",abs=True)
pypath    = dirpath('soft/anaconda/bin/')
python    = exepath(pypath+'python',abs=True)
phonopy   = exepath(pypath+'phonopy',abs=True)
phono3py  = exepath(pypath+'phono3py',abs=True)
phonts    = exepath('soft/PhonTS-1.1.4/src/PhonTS')
vasp      = exepath('vasp')
vasppot   = dirpath('../zhangyueyu/psudopotential/5.3_LDA')
shengbte  = exepath('tcscripts/ShengBTE/ShengBTE')
thirdorder= python+exepath('soft/thirdorder/python/thirdorder_vasp.py')
libs	  =['/opt/intel/mkl/10.0.013/lib/em64t','/opt/intel/mpi/openmpi/1.6.3/icc.ifort/lib']