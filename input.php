<?php
require_once("discript.php");
#require_once("postMini.php");
if(!$langevin)$langevin=0;
if(!$nvt)$nvt=0;
if(!$jprofile)$jprofile=0;
if(!$computeTc)$computeTc=0;
if(!$fourierTc)$fourierTc=0;
if(!$gstart)$gstart=0;
if(!$jcf)$jcf=0;
if(!$dumpxyz)$dumpxyz=0;
if(!$dumpv)$dumpv=0;#$potential='b';$masses="b";
if(!$hdeta)$hdeta=0.0
$kb=$boltz[$units];
$nktv=$nktv2p[$units];
$bb="python $home/input.py $units $xp $yp $zp $dumpRate $timestep $method $kb $nktv '$masses' '$potential' $T $seed $dtime $equTime $langevin $nvt $aveRate $deta $jprofile $dumpRate $corRate $computeTc  $fourierTc $tcfactor  $gstart $jcf  $nswap $excRate  $excNum $swapEnergyRate $dumpxyz $dumpv $runTime $upP $wfix $nstat '$enforceThick' '$thick' $Thi $Tlo '$hdeta'";
passthru($bb);

