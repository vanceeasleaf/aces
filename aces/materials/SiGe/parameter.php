<?php
if(!$units)$units="metal";
//$species="CNSi";//or C
//***********general***********
if(!$method)$method="inject";
$enforceThick=0;
$yp=1;#yPeriodic
$zp=1;
$useMini=1;
$dumpxyz=1;
$upP=3;
//in metal default
$srcUnit="metal";
$usinglat=1;
$xlen=100*Unit($srcUnit,"l");
$latx=288;$laty=4;$latz=4;
$ylen=20*Unit($srcUnit,"l");
$thick=3.35*Unit($srcUnit,"l");
$deta=2.4*Unit($srcUnit,"l");
$T=500*Unit("metal","T");//K for all units
$Thi=$T+10*Unit("metal","T");
$Tlo=$T - 10*Unit("metal","T");

$timestep=.55e-3*Unit($srcUnit,"t");
$equTime=100000;
$dumpRate=10000;
$aveRate=100000;
$excRate=1;
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
$bond=5.430*Unit($srcUnit,"l");
$ratio=0.1;
?>