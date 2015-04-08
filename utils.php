<?php
/* 程序中用到的所有文件名*/
$fileTempProfile="$projHome/tempProfile.txt";
$fileTempAve="$projHome/tempAve.txt";
$fileXRange="$projHome/minimize/minimize.xyz";
$fileNvtWork="$projHome/nvtWork.txt";
$fileSwap="$projHome/swapEnergy.txt";
$fileBoxRange="$projHome/minimize/range";
$fileRestartMini="$projHome/minimize/restart.minimize";
$fileDump="$projHome/dump.lammpstrj";
$fileDumpVelocity="$projHome/dump.velocity";
$fileTherCon="$projHome/thermal_conductivity.txt";
$fileTempDiff="$projHome/delta_temp.txt";
$fileKappa="$projHome/kappa.txt";
$fileFlux="$projHome/flux.txt";
$fileResult="$projHome/result.txt";
$fileScanResult="$projHome/scan.txt";
$fileJCF="$projHome/jcf.txt";
$fileJProfile="$projHome/jprofile.txt";
$dtime=$timestep*100;
$excNum=$aveRate/$excRate;
$corNum=$aveRate/$corRate;
$swapEnergyRate=$swapEnergy/($excRate*$timestep);
$xp=1;#xPeriodic
if($method=="nvt")$xp=0;

require_once("$home/materials/".$species."/disp.php");
$home=dirname(__FILE__);	
require_once("$home/funcs.php");

/* 将热导率转化为国际单位的因子*/
$tcfactor=($boltz["si"]*$kelvin["si"])/($boltz[$units]*$kelvin[$units])*$femtosecond[$units]/$femtosecond["si"]*$angstrom[$units]/$angstrom["si"]*$kelvin[$units]/$kelvin["si"];//the unit is J/s/m/K

    
?>
