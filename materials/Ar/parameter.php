<?php
if(!$units)$units="lj";
//$species="CNSi";//or C
//***********general***********
if(!$method)$method="muller";
$enforceThick=0;
$yp=1;#yPeriodic
$zp=1;
$useMini=0;
$dumpxyz=1;
$upP=3;
//in metal default
$usinglat=1;
$srcUnit="lj";
$xlen=2*.844*8.000000*Unit($srcUnit,"l");
$ylen=2*.844*8.000000*Unit($srcUnit,"l");
$thick=2*.844*8*Unit($srcUnit,"l");
$deta=2*.844*.25*Unit($srcUnit,"l");
$latx=8;$laty=8;$latz=8;
$T=.71*Unit("lj","T");//K for all units
$Thi=$T+1.0*Unit("lj","T");
$Tlo=$T - 1.0*Unit("lj","T");

$timestep=0.000466*Unit($srcUnit,"t");
$equTime=100000;
$dumpRate=10000;
$aveRate=100000;
$excRate=10;
$corRate=2;
$runTime=1000000;
$seed=458127641;
//***********for nvt***********
$jcf=1;
$nstat=1;//heat bath width;unit=$deta
$nswap=$nstat;
$wfix=3;//fix width
$Nbins=500;
$begin=1;
$cutoff=2.8*Unit("lj","l");
$epsilon1=1*Unit("lj","E");
$sigma1=1*Unit("lj","l");
$latticeConstant=.844*Unit($srcUnit,"l")*($units!="lj"?2:1);
$m1=1*Unit("lj","M");
?>