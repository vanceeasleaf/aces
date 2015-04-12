<?php
require_once("funcs.php");
function genPbss($path,$disp,$queue,$nodes,$procs,$start,$ucores){
$len=floor($nodes*$procs/$ucores);
$real_cores=$len*$ucores;
$end=$start+$len-1;
	require("config.php");
	if(!is_dir("$path/pbs/minimize"))mkdir("$path/pbs/minimize",0777,true);
$pbsfile="$path/pbs/lammps${start}-${end}.pbs";
$pbs=fopen($pbsfile,"w");
if(!is_file($pbsfile))die("create file  $pbsfile failed!\n");
$mfile="$path/pbs/minimize/in.${start}-${end}";
$input_m=fopen($mfile,"w");
$ifile="$path/pbs/in.${start}-${end}";
$input=fopen($ifile,"w");$b=$len;
for($i=$start;$i<=$end;$i++){
	$p=$i-$start+1;
	fprintf($input_m,"
	partition yes $p log $path/$i/minimize/log.lammps
	partition yes $p jump $path/$i/minimize/input
	");
	fprintf($input,"
	partition yes $p  log $path/$i/log.lammps
	partition yes $p 	jump $path/$i/input
	");
}
#define MPI PATH
	$home=dirname(__FILE__);	
	$php="$PHP_HOME/php";
fprintf($pbs,
"#!/bin/bash -x
#PBS -l nodes=$nodes:ppn=$procs
#PBS -l walltime=240:00:00
#PBS -j oe
#PBS -q $queue
#PBS -N ${disp}${start}-${end}


# Setup the OpenMPI topology
n_proc=\$(cat \$PBS_NODEFILE | wc -l)

contexts=`~/bin/get_psm_sharedcontexts_max.sh`
 if [ '$?' = '0' ] ; then
  export PSM_SHAREDCONTEXTS_MAX=$contexts
 fi
");
for($i=$start;$i<=$end;$i++){
	fprintf($pbs,
	"cd  $path/$i/minimize
	$php $home/minimize/input.php \"$path/$i/qloop.php\" \"$path/$i/species.php\">input  
	");
}
fprintf($pbs,
"	$OMPI_HOME/bin/mpirun  -machinefile \$PBS_NODEFILE -np $real_cores  $APP_PATH -plog $path/pbs/minimize/log.${start}-${end} -partition ${len}x${ucores} -pscreen $path/pbs/minimize/screen.${start}-${end} -i $mfile
	");
for($i=$start;$i<=$end;$i++){
	fprintf($pbs,
	"cd  $path/$i
	$php $home/input.php  \"$path/$i/qloop.php\" \"$path/$i/species.php\">input  
	");
}
fprintf($pbs,
"$OMPI_HOME/bin/mpirun  --mca btl openib,self  -machinefile \$PBS_NODEFILE -np $real_cores   $APP_PATH -partition ${len}x${ucores} -plog $path/pbs/log.${start}-${end} -pscreen $path/pbs/screen.${start}-${end} -i $ifile
");
fprintf($pbs,
"exit 0");
}
?>