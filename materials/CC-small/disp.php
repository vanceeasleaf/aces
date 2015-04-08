<?php
$potential="
pair_style        tersoff
pair_coeff      * * $home/potentials/BNC.tersoff  C C
";
$dump="
dump_modify dump1 element C D
";
$home=dirname(__FILE__);
require_once("$home/structure.php");
function structure(){
	global $projHome;
	$home=dirname(__FILE__);
	require_once("$home/../../config.php");
preFile();
	echo"
read_data $projHome/minimize/structure
";
global $mass;
if($mass)echo "
	mass 2 $mass
";
}
?>