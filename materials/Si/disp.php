<?php

$potential="
pair_style        sw
pair_coeff      * * $home/potentials/Si.sw  Si
";
$dump="
dump_modify dump1 element Si
";
function structure(){
	$home=dirname(__FILE__);
	require_once("$home/../../config.php");
	$php="$PHP_HOME/php";
 shell_exec("$php $home/structure.php");
	echo"
read_data structure
";
}
?>