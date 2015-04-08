<?php
if(!$xp||!yp||!zp)exit("Material Ar only support totally periodic environment!\n");
$potential="
pair_style        lj/cut $cutoff
pair_coeff        1  1   $epsilon1 $sigma1    #  LJ parameters for Ar-Ar 
";
$dump="
dump_modify dump1 element Ar
";

function structure(){
	$home=dirname(__FILE__);
	require_once("$home/../../config.php");
	$php="$PHP_HOME/php";
 echo shell_exec("$php $home/structure.php");
}
$masses="
mass 1 $m1
";
?>