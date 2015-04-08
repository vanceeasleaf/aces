<?php
	$home=dirname(__FILE__);	
require_once("$home/../../discript.php");

?>
units            <?echo $units?> 
atom_style atomic
boundary p p p
dimension 3
read_data rand_structure
<?if($ratio===0){
echo "
	set group all type 1
";
}else if($ratio===1){
	echo "
	set group all type 2
";
}
?>
<?	echo $potential;?> 
<?	echo $masses;?> 
thermo_style custom step pe etotal
thermo 1
dump 2 all xyz  1 metropolis.xyz
dump_modify 2 element C N
min_style metropolis
minimize 1e-6 1e-12 1000000 1000000
write_data structure