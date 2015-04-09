<?php
/**
 * MD生成的文件的后处理
 * @author zhouy
 */
 
require_once("discript.php");
require_once("postMini.php");

if(!$conti)$conti=0;
if(!$fourierTc)$fourierTc=0;
if(!$computeTc)$computeTc=0;
$kb=$boltz[$units];
passthru("cd $projHome;python $home/profile.py $method  $begin $timestep $S  $conti $lz $excRate $swapEnergyRate $upP $deta $tcfactor $zfactor $lx $fourierTc $computeTc $corRate $kb $ly $T");



?>