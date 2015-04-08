<?php
	$fp=fopen("metropolis.structure","r");
for($i=0;$i<3;$i++)
fgets($fp);
list($natom)=fscanf($fp,"%d");
fgets($fp);
list($xlo,$xhi)=fscanf($fp,"%f%f");
list($ylo,$yhi)=fscanf($fp,"%f%f");
list($zlo,$zhi)=fscanf($fp,"%f%f");
fgets($fp);
for($i=0;$i<$natom;$i++){
	list($id,$type[$i],$xs[$i],$ys[$i],$zs[$i])=fscanf($fp,"%d%d%f%f%f");
}
for($i=0;$i<3;$i++)
fgets($fp);
list($natom)=fscanf($fp,"%d");
fgets($fp);
list($xlo,$xhi)=fscanf($fp,"%f%f");
list($ylo,$yhi)=fscanf($fp,"%f%f");
list($zlo,$zhi)=fscanf($fp,"%f%f");
fgets($fp);
for($i=0;$i<$natom;$i++){
	list($id,$type[$i],$xs[$i],$ys[$i],$zs[$i])=fscanf($fp,"%d%d%f%f%f");
}
$fout=fopen("structure","w");
fprintf($fout,"# Graphene cell 
       $natom  atoms
         2  atom types

 $xlo    $xhi  xlo xhi
$ylo     $yhi  ylo yhi
$zlo     $zhi  zlo zhi

Masses

1 12.0112
2 14.0067

Atoms

");
for($i=0;$i<$natom;$i++){
	fprintf($fout,"%d %d     %f     %f     %f\n",$i+1,$type[$i],$xs[$i]*($xhi-$xlo),$ys[$i]*($yhi-$ylo),$zs[$i]*($zhi-$zlo));
}