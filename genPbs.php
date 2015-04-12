<?php
require_once("funcs.php");
function genPbs($path,$disp,$queue,$nodes,$procs){

	require("config.php");
$pbs=fopen("$path/lammps.pbs","w");
#define MPI PATH
	$home=dirname(__FILE__);	
	$php="$PHP_HOME/php";
$LOGOUT="log.out";
fprintf($pbs,
"#!/bin/bash -x
#PBS -l nodes=$nodes:ppn=$procs
#PBS -l walltime=240:00:00
#PBS -j oe
#PBS -q $queue
#PBS -N $disp


# Setup the OpenMPI topology
n_proc=\$(cat \$PBS_NODEFILE | wc -l)

contexts=`~/bin/get_psm_sharedcontexts_max.sh`
 if [ '$?' = '0' ] ; then
  export PSM_SHAREDCONTEXTS_MAX=$contexts
 fi
cd  $path/minimize
$php $home/minimize/input.php \"$path/qloop.php\" \"$path/species.php\">input   
$OMPI_HOME/bin/mpirun  -machinefile \$PBS_NODEFILE -np \$n_proc $APP_PATH <input &>$path/minimize/$LOGOUT
cd  $path
$php $home/input.php \"$path/qloop.php\" \"$path/species.php\">input
$OMPI_HOME/bin/mpirun  --mca btl openib,self  -machinefile \$PBS_NODEFILE -np \$n_proc $APP_PATH <input &>$path/$LOGOUT


exit 0");
}
function genSh($path,$disp,$procs){
		require("config.php");
		$pbs=fopen("$path/run.sh","w");
			$home=dirname(__FILE__);	
	$php="$PHP_HOME/php";
$LOGOUT="log.out";
fprintf($pbs,
"#!/bin/bash -x
#$disp
cd  $path/minimize
$php $home/minimize/input.php \"$path/qloop.php\" \"$path/species.php\">input   
mpirun   -np $procs $APP_PATH <input >$path/minimize/$LOGOUT 2>/dev/null
cd  $path
$php $home/input.php \"$path/qloop.php\" \"$path/species.php\">input
mpirun   -np $procs $APP_PATH <input >$path/$LOGOUT 2>/dev/null
exit 0");
}
?>