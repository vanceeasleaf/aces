<?php
$potential="
pair_style        tersoff
pair_coeff      * * $home/potentials/BNC.tersoff  C N
";
$dump="
dump_modify dump1 element C N
";
function structure(){
	global $projHome;
	$home=dirname(__FILE__);
	require_once("$home/../../config.php");
preFile();
	echo"
read_data $projHome/minimize/structure
		mass 1 12.01
	mass 2 14.00
";
}

/* C2N hollow 2D structure generator.*/ 

function preFile(){
	global $laty;
	global $latx;
	global $sideLen;/* the number of Nitrogen on each side*/
	if($latx%2==1)$latx+=1;
      global $bond;

 $home=dirname(__FILE__);	
passthru("python $home/structure.py $latx $laty $bond $sideLen");

}
?>