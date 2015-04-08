<?php
if(!$units)$units="metal";
//$species="CNSi";//or C
//***********general***********
if(!$method)$method="nvt";
$enforceThick=1;
$yp=0;#yPeriodic
$zp=1;
$useMini=1;
$dumpxyz=1;
$upP=3;
//in metal default
$srcUnit="metal";
$usinglat=0;
$xlen=100*Unit($srcUnit,"l");
$latx=228;$laty=4;$latz=4;
$ylen=20*Unit($srcUnit,"l");
$thick=3.35*Unit($srcUnit,"l");
$deta=2.5*Unit($srcUnit,"l");
$T=300*Unit("metal","T");//K for all units
$Thi=$T+10*Unit("metal","T");
$Tlo=$T - 10*Unit("metal","T");

$timestep=.182e-3*Unit($srcUnit,"t");
$equTime=100000;
$dumpRate=100000;
$aveRate=100000;
$excRate=500;
$corRate=2;
$runTime=10000000;
$seed=458127641;
//***********for nvt***********

$nstat=1;//heat bath width;unit=$deta
$nswap=$nstat;
$swapEnergy=5e-4*Unit("metal","E");
$wfix=2;//fix width
$Nbins=500;
$bond=1.42*Unit($srcUnit,"l");
$begin=1;
$cell="cell.in";
	
$latysmall=3; /* 窄边包括2*$latysmall+1排C原子，$latysmall=3时，为7 AGNR*/
$latExtend=2;/* 相对于窄边，长边的y方向两头都延伸$latExtend个六边形，$latExtend=2时，长边为[($latysamll-1)+$latExtend]*2+1=13 AGNR*/
$latxsmall=7;//奇数
$latxbig=12;/* 最简单模型，中间是长边，两头是等长的短边, $latxsmall表示短边由几条*折线*组成，偶数*/

?>