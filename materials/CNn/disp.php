<?php
if($ratio==0)$NN="";else $NN="N";
$potential="
pair_style        tersoff
pair_coeff      * * $home/potentials/BNC.tersoff  C $NN
";
$dump="
dump_modify dump1 element C $NN
";
function structure(){
	$home=dirname(__FILE__);
	require_once("$home/../../config.php");
	$php="$PHP_HOME/php";
echo shell_exec("$php $home/structure.php");
	echo"
read_data structure
";
}
?>