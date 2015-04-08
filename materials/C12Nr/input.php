<?php
	$home=dirname(__FILE__);	
require_once("$home/../../discript.php");

?>
units            <?echo $units?> 
atom_style atomic
boundary p p p
dimension 3
read_data rand_structure
<?	echo $potential;?> 
<?	echo $masses;?> 
thermo_style custom step pe etotal
thermo 1
dump 1 all atom 1 metropolis.lammpstrj
dump 2 all xyz  1 metropolis.xyz
dump_modify 2 element C N
dump 3 all atom 10000000 metropolis.structure
min_style metropolis
minimize 1e-12 1e-12 1000000 1000000