<?php

if(!$enforceThick)exit("The problem is involving a sheet ,please set enforceThick true!\n");
$potential="
pair_style        tersoff
pair_coeff      * * $home/potentials/BNC.tersoff  C N  C
";
$dump="
dump_modify dump1 element C N  Si
";
function structure(){
		$home=dirname(__FILE__);	
require_once("$home/../../config.php");
	$php="$PHP_HOME/php";
echo shell_exec("$php $home/structure.php");
echo "
read_data structure
";
}
$masses="
mass 1 $m1
mass 2 $m2
mass 3 $m3
";
?>