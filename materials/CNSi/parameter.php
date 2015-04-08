<?php
if(!$units)$units="metal";
//$species="CNSi";//or C
//***********general***********
if(!$method)$method="muller";
$enforceThick=1;
$yp=1;#yPeriodic
$zp=1;
$useMini=1;
$dumpxyz=1;
$upP=3;
//in metal default
$srcUnit="metal";
$xlen=100*Unit($srcUnit,"l");
$ylen=20*Unit($srcUnit,"l");
$thick=3.35*Unit($srcUnit,"l");
$deta=3*Unit($srcUnit,"l");
$T=300*Unit("metal","T");//K for all units
$Thi=$T+10*Unit("metal","T");
$Tlo=$T - 10*Unit("metal","T");

$timestep=.182e-3*Unit($srcUnit,"t");
$equTime=100000;
$dumpRate=10000;
$aveRate=100000;
$excRate=500;
$corRate=2;
$runTime=10000000;
$seed=458127641;
//***********for nvt***********

$nstat=3;//heat bath width;unit=$deta
$nswap=$nstat;
$wfix=3;//fix width
$Nbins=500;

$log=1;
$xyz=1;
$xZigzag=1;
$rmass=0;//ratio of strange C
$cn=0;
$bond=1.42*Unit("metal","l");
$m1=12.01*Unit("metal","M");
$m2=14*Unit("metal","M");
$m3=32*Unit("metal","M");
?>