<?php

if(!$enforceThick)exit("The problem is involving a sheet ,please set enforceThick true!\n");
$potential="
pair_style        tersoff
pair_coeff      * * $home/potentials/SiC_1989.tersoff  C C C
";
$dump="
dump_modify dump1 element C  C C
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
mass 1 12.01
mass 2 12.01
mass 3 12.01
";
?>
