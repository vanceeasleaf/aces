<?php
$potential="
pair_style        tersoff
pair_coeff      * * $home/potentials/BNC.tersoff  B N
";
$dump="
dump_modify dump1 element B N
";
$home=dirname(__FILE__);
require_once("$home/structure.php");
function structure(){
	$home=dirname(__FILE__);
	require_once("$home/../../config.php");
preFile();
	echo"
read_data structure
";
}
?>