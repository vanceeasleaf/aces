#!/bin/bash
cd $PBS_O_WORKDIR
module load openmpi-psm-gcc

mpirun -np `cat $PBS_NODEFILE | wc -l` EXECPATH/EXECUTABLE < RUNPATH/LMP_TEMP 

