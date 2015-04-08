<?php
echo testEnergy(1,1,0,1.42)."\n";
echo testEnergy(1,2,0,1.42)."\n";
echo testEnergy(2,2,0,1.42)."\n";
function testEnergy($type1,$type2,$x,$y){
$file=fopen("bond","w");
fprintf($file,
"# bond
       2  atoms
         2  atom types

 0.    5  xlo xhi
 0.    5  ylo yhi
 0.    5  zlo zhi

Masses

1 12.0112
2 14.0067

Atoms

1 $type1    0    0    0
2 $type2    $x     $y     0
");
$energy=shell_exec("~/home1/zhouy/lmp_ubuntu <in.bond >err 2>err;tail -22 log.lammps|head -1");
return trim($energy);
}