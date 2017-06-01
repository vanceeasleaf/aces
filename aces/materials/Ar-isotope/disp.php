<?php
//if(!$xp||!yp||!zp)exit("Material Ar only support totally periodic environment!\n");
$potential="
pair_style        lj/cut $cutoff
pair_coeff        *  *   $epsilon1 $sigma1    #  LJ parameters for Ar-Ar 
";
$dump="
dump_modify dump1 element Ar Ap
";
	$home=dirname(__FILE__);
function structure(){
		global $projHome,$ratio,$usinglat,$latx,$laty,$latz,$xlen,$ylen,$thick,$bond,$latticeConstant,$m2;
$home=dirname(__FILE__);
	require_once("$home/../../config.php");
	$php="$PHP_HOME/php";
 if(!$usinglat){
$latx=floor($xlen/(2*$bond));
$laty=floor($ylen/(2*$bond));
$latz=floor($thick/(2*$bond));
	}
	echo "
lattice            fcc    $latticeConstant
region            box  block 0  $latx  0 $laty  0 $latz  units lattice
create_box     2  box
create_atoms     1  box
set group all  type/fraction 2 $ratio 12393
";
}
$masses="
mass 1 $m1
mass 2 $m2
";
?>