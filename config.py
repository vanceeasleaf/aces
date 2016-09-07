root="/home/xggong/home1/zhouy/"
def exepath(a,abs=False):
	if abs:
		return ' '+a+' '
	else:
		return ' '+dirpath(a)+' '
def dirpath(a):
	return root+a.lstrip('/')

php       = exepath("php/bin/php")
phana     = exepath("soft/lammps-25Sep15/tools/phonon/phana")
espresso  = dirpath('soft/QE+EPW/espresso-4.0.3/')
epw		  = exepath(espresso+"EPW/bin/epw.x",abs=True)
ph        =	exepath(espresso+"/bin/ph.x",abs=True)
pw        =	exepath(espresso+"/bin/pw.x",abs=True)
boltztrap = exepath("soft/boltztrap-1.2.5/src/BoltzTraP")
x_trans = exepath("soft/boltztrap-1.2.5/src/x_trans")
boltztrap_vasp = exepath("soft/boltztrap-1.2.5/zhangwenqing/BoltzTrap_vasp/BoltzTrap_vasp")
x_trans_vasp = exepath("soft/boltztrap-1.2.5/zhangwenqing/BoltzTrap_vasp/x_trans")
phanalib  = dirpath("soft/phona")
lammps    = exepath("soft/lammps-25Sep15/src/lmp_mpi")
lammpspot = dirpath("soft/lammps-25Sep15/potentials")
mpirun    = exepath("/opt/intel/mpi/openmpi/1.6.3/icc.ifort/bin/mpirun  -np ",abs=True)
pypath    = dirpath('soft/anaconda/bin/')
python    = exepath(pypath+'python',abs=True)
phonopy   = exepath(pypath+'phonopy',abs=True)
phono3py  = exepath(pypath+'phono3py',abs=True)
phonts    = exepath('soft/PhonTS-1.1.4/src/PhonTS')
vasp      = exepath('vasp')
vasppot   = dirpath('../zhangyueyu/psudopotential')
shengbte  = exepath('tcscripts/ShengBTE/ShengBTE')
thirdorder= python+exepath('soft/thirdorder/thirdorder_vasp.py')
alamode= dirpath('soft/alamode-develop/')
alm=exepath(alamode+'alm/alm',abs=True)
anphon=exepath(alamode+'anphon/anphon',abs=True)
almdisp=python+exepath(alamode+'tools/displace.py --VASP=POSCAR-supercell ',abs=True)
#thirdorder= python+exepath('soft/thirdorder0.9/python/thirdorder_vasp.py')
libs	  =['/opt/intel/mkl/10.0.013/lib/em64t','/opt/intel/mpi/openmpi/1.6.3/icc.ifort/lib']