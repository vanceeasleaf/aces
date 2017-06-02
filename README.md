[![Build Status](https://travis-ci.org/vanceeasleaf/aces.svg?branch=master)](https://travis-ci.org/vanceeasleaf/aces)
[![PyPi version](https://img.shields.io/pypi/v/aces.svg)](https://pypi.python.org/pypi/aces)
[![PyPi license](https://img.shields.io/pypi/l/aces.svg)](https://pypi.python.org/pypi/aces)
[![PyPi downloads](https://img.shields.io/pypi/downloads/aces.svg)](https://pypi.python.org/pypi/aces)
[![Python version](https://img.shields.io/pypi/pyversions/aces.svg)](https://pypi.python.org/pypi/aces)

# Automatical Computational Experiment System(ACES)

A python framework for **computational physics** numerical experiments.

> @Authour Yang Zhou

> @Mail  y_zhou13@fudan.edu.cn

> @Fudan University Computational Condensed Matter Group(CCMG)

Docs at [http://vanceeasleaf.github.io/aces/](http://vanceeasleaf.github.io/aces/)

# Summary

## What

- ACES is a wrapper for many computational codes for **thermal conductivity** with similar workflows.

- More generally, it provides a framework to automate the computing process including pre-process and post-process.

## Why

- First principle calculation or molecular dynamics computation code usually needs structure files such as POSCAR of **VASP**,  data file of **LAMMPS**, struct file of **Quantum Espresso** and  input files such as INCAR, KPOINTS of **VASP**, CONTROL file of **ShengBTE**. The information of these files are similar but the format needs carefully treatment,for instance, the length and energy units. ACES is pretended to unify the input and out file format of all the materials computational codes to let user easily carry up a new engine.

- A computational experiment project is always  chaotic because of confused folder structure and name. ACES let you manage your project in a unify standard which could be version controlled and is easily to track.

- The software engineer method is applyed to computational experiment to extract reuseable materials and devices which let you contribute your module to the public or enjoy the materials database from  other scientists.

## How

- Runners are defined to wrap a certain computational code and are design to be assemble of micro functional modules which could be used sequenly to complete a certain job or be used though command-line to controll the job flexiblely.

- ACES uses json to manage the options of a project . A json file consists of  materials to be computed,  runner  used to compute it,  computing parameters such as temperature and sample length, and how many cores are intented to be used.

- A set of thermal conductivity calculation algorithm including NEMD, Muller-Plathe, Green-Kubo, NEGF, NMA, BTE are realized.  

- Input files and structure files for LAMMPS, VASP, ShengBTE, alamode etc. are ready to be generated.

- PBS are used to manage your tasks. 

- Automatically analysis including data formatting and graph generation are carried out.

- Very easily for you to extend ACES if you have some new codes to integrate. For example, you can write a wrapper for Quantum Expresso easily to use it as an engine to structure optimization, band structure calculation or phonon calculation to have the same work flow with VASP.


# Usage

## Requirements

We recommand you to install Anaconda-2.4.sh to prepare a python2.7 environment.

requrements list:

- ase==3.13.0
- numpy==1.12.1
- Cython==0.25.2
- Image==1.5.5
- Numeric==24.2
- PyOpenGL==3.1.1a1
- Pmw==2.0.1
- h5py==2.7.0
- lxml==3.7.3
- matplotlib==2.0.2
- mpi4py==2.0.0
- np==0.2.0
- pandas==0.20.1
- pexpect==4.2.1
- pyspglib==1.8.3.1
- scanf==1.4.1
- scipy==0.19.0
- scikit_learn==0.18.1
- atlas==0.27.0
- vapory==0.1.01
- PyYAML==3.12


## Install 

### Online install

You can easily install ACES by 

``` bash 
pip install aces
```

### Manually install

After downloading you can run 

``` bash
unzip aces-master.zip 
cd aces-master
pip install -r requirements.txt
```
to install the dependecies.


If you are offline, you can download the reqirements by 

``` bash 
pip install -r requirements.txt --find-links=.
```
and all the dependencies will be downloaded to the current directory. 

After copy them to your target computer you can run 
``` bash 
pip install -r requirements.txt --find-links=.
```
to install the dependecies.

## Basic Usage

run 
``` bash
python setup.py install
```
to install it to the %python%/site-packages/aces folder.

run 
``` bash
python -c 'import aces'
``` 
to check if the installation is successful.

You have to change the path in config.py to tell ACES where is LAMMPS, ShengBTE,Phonopy and so on.

### Example 
Let's calculate the phonon dispertion, thermal conductivity and so on of black phosphorous using ACES with VASP engine and Phonopy, ShengBTE. 

in your workspace 

``` bash 
mkdir BP.project 
cd BP.project 
touch sub.py
```

input `sub.py`

``` python
from aces import Aces
#the origin BP structure is optimized and we use it directly
class sub(Aces):
	def submit(self):
		opt=dict(
			units="metal",
			species="BP",
			method="nvt",
			nodes=1,
			procs=12,
			queue="q1.4",
			runTime=10000000
			,runner="shengbte"
		)
		app=dict(shengcut=-4,
		kpoints=[61,61,1],
		engine='vasp',
		supercell=[4,4,1],
		ekpoints=[1,1,1])
		self.commit(opt,app);
if __name__=='__main__':
	sub().run()
```

Then run 
``` bash 
python sub.py
``` 
The only thing you have to do is waiting for the result.

Detailed document at at [http://vanceeasleaf.github.io/aces/](http://vanceeasleaf.github.io/aces/)
