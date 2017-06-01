<?php
$potential="
pair_style        tersoff
pair_coeff      * * $home/potentials/SiCGe.tersoff  Si(D) Ge
";
$dump="
dump_modify dump1 element Si Ge
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
}
?>