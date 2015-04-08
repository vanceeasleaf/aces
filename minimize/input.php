<?php
	$home=dirname(__FILE__);	
require_once("$home/../discript.php");

?>
units            <?echo $units?> 
atom_style atomic
boundary p p p
dimension 3
<? structure();?> 
<?	echo $potential;?> 
timestep          <?echo $timestep?> 
<?	echo $masses;?> 
thermo_style custom step pe etotal
thermo <?echo $dumpRate?> 
<?if($write_structure){
	echo "write_data structure
	dump dumpc all xyz 1 CN.xyz
	run 0 ";
}?>
<?if($metropolis){
echo"
min_style metropolis
minimize 1e-12 1e-12 1000000 1000000
";
}
?> 
<?if($useMini){
echo"
fix 1 all box/relax x 0.0 y 0.0 nreset 1
#dump 1 all custom 10 dump.minimize.* type id x y z
#dump Graphene all xyz 1000000 graphene.*.xyz
#dump_modify Graphene element C N C

min_style cg
minimize 1e-12 1e-12 1000000 1000000
";
}
?> 
write_restart <?echo $fileRestartMini?> 
dump dump1 all xyz 1 <?echo $fileXRange?> 
<?	echo $dump;?> 
dump kaka all atom 1 <?echo $fileBoxRange?> 
run 0
#fix frelax all nve
#fix controltemp all temp/rescale 100 10.0 10.0 10.0 1.0
#run 1000
#unfix controltemp
#unfix frelax
#undump Graphene
