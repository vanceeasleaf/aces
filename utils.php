<?php
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
