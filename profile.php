<?php
/**
 * MD生成的文件的后处理
 * @author zhouy
 */
 
require_once("discript.php");

if(!$conti)$conti=0;
if(!$fourierTc)$fourierTc=0;
if(!$computeTc)$computeTc=0;
$kb=$boltz[$units];
//echo "cd $projHome;/home/xggong/home1/zhouy/soft/anaconda/bin/python $home/profile.py $method  $begin $timestep   $conti  $excRate $swapEnergyRate $upP $deta $tcfactor  $fourierTc $computeTc $corRate $kb  $T $xp $yp $zp $enforceThick $thick";
passthru("cd $projHome;/home/xggong/home1/zhouy/soft/anaconda/bin/python $home/profile.py $method  $begin $timestep   $conti  $excRate $swapEnergyRate $upP $deta $tcfactor  $fourierTc $computeTc $corRate $kb  $T $xp $yp $zp $enforceThick $thick");



?>