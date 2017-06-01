<?php
$units="metal";
//$species="CNSi";//or C
//***********general***********
$method="muller";
$enforceThick=0;
$yp=1;#yPeriodic
$zp=1;
$useMini=1;
$dumpxyz=1;
$upP=3;
//in metal default
$srcUnit="metal";
$usinglat=1;
$xlen=40*Unit($srcUnit,"l");
$latx=5;$laty=5;$latz=2;
$ylen=40*Unit($srcUnit,"l");
$thick=40*Unit($srcUnit,"l");
$deta=1.4*Unit($srcUnit,"l");
$T=300*Unit("metal","T");//K for all units
$Thi=$T+10*Unit("metal","T");
$Tlo=$T - 10*Unit("metal","T");

$timestep=.182e-3*Unit($srcUnit,"t");
$equTime=100000;
$dumpRate=10000;
$aveRate=100000;
$excRate=500;
$corRate=2;
$runTime=2000000;
$seed=458127641;
//***********for nvt***********

$nstat=3;//heat bath width;unit=$deta
$nswap=$nstat;
$swapEnergy=5e-4*Unit("metal","E");
$wfix=3;//fix width
$Nbins=500;

$begin=1;
?>