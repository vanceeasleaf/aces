<?php
/* C2N triangle structure generator.*/ 

/**
 * rotate around z axis by $phi
 * @author zhouy
 */
function rotate($phi,$x,$y){
	$x1=cos($phi)*$x-sin($phi)*$y;
	$y1=sin($phi)*$x+cos($phi)*$y;
	return array($x1,$y1);
}

/**
 * rotate a group of atoms
 * @author zhouy
 * @input $atoms a 2D array,$phi the degree
 */
function rotateAtom($phi,$atoms){
	$Natom=count($atoms);
	$newAtoms=array();
	for($i=0;$i<$Natom;$i++){
		$a=$atoms[$i];
		$b=rotate($phi,$a[0],$a[1]);
		array_push($b,$a[2]);
		array_push($newAtoms,$b);
	}
	return $newAtoms;
}

/**
 * move a group of atoms
 * @author zhouy
 */
function moveAtom($x,$y,$atoms){
	$Natom=count($atoms);
	$newAtoms=array();
	for($i=0;$i<$Natom;$i++){
		$a=$atoms[$i];
		$b=array($a[0]+$x,$a[1]+$y,$a[2]);
		array_push($newAtoms,$b);
	}
	return $newAtoms;
}
/**
 * generate a 9-atom unit cell
 * @author zhouy
 */
function getUnitCell(){
	
	/* minimal generator*/
	$atoms0=array(
		array(1/2*sqrt(3),0.5,'C'),
		array(1/2*sqrt(3),-0.5,'C'),
		array(1/2*sqrt(3),1.5,'N')
		);
	$unitCell=array();
	for($i=0;$i<3;$i++){
		$newAtoms=rotateAtom(2*3.1415926/3*$i,$atoms0);
		$unitCell=array_merge($unitCell,$newAtoms);
	}
	global $cen;
	if($cen)array_push($unitCell,array(0,0,'N'));
	return $unitCell;
}

/** 
 * the triangle structure doesn't have a rectangular cell, a quasi rectangular cell is therefore built
 * @author zhouy
 * @input $row=y ,$col=x,$type=the most left is above if 0
 */
function quasiCell($row,$col,$type=0){
	$atoms0=getUnitCell();
	$unitCell=array();
	for($i=0;$i<$col;$i++){
		for($j=0;$j<$row;$j++){
			$x=$i*6*sqrt(3)/2;
			$y=$j*6;
			$newAtoms=moveAtom($x,$y,$atoms0);
			$unitCell=array_merge($unitCell,$newAtoms);
		}
	}
	return $unitCell;
}
function preFile(){
	global $laty;
	global $latx;
      global $bond;
	$unitCell=quasiCell($laty,$latx);
    	$Natom=count($unitCell);
    	$x=array();
    	$y=array();
    	$la=array();
    	/* scale and translate the coordinates*/
    	for($i=0;$i<$Natom;$i++){
    		$x[$i]=$bond*$unitCell[$i][0];
    		$y[$i]=$bond*$unitCell[$i][1];
    		$la[$i]=$unitCell[$i][2];
    	}
    	$minx=$miny=10000;
    	$maxx=$maxy=-10000;
    	for($i=0;$i<$Natom;$i++){
    		$minx=min($minx,$x[$i]);
    		$miny=min($miny,$y[$i]);
    		$maxx=max($maxx,$x[$i]);
    		$maxy=max($maxy,$y[$i]);
    	}
    	$xs=array();
    	$ys=array();
    	$length=$latx*6*$bond*sqrt(3)/2;
    	$width=$laty*6*$bond;
   	for($i=0;$i<$Natom;$i++){
    		$xs[$i]=$x[$i]/$length;
    		$ys[$i]=$y[$i]/$width;
    	}
    	
    	/* sort the atoms*/
    	$cs=array();
    	$ns=array();
    	for($i=0;$i<$Natom;$i++){
    		if($la[$i]=='C')array_push($cs,$i);
    		if($la[$i]=='N')array_push($ns,$i); 			
    	}
	$nC=count($cs);
	$nN=count($ns);    	
/* generate cell.in*/
$file=fopen("cell.in","w");
fprintf($file,"1.0
$length 0 0
0 $width 0 
0 0 100.0
$nC $nN
");
for($j=0;$j<$nC;$j++){
	$i=$cs[$j];
	fprintf($file,"%f\t%f\t0.5\n",$xs[$i],$ys[$i]);
}
for($j=0;$j<$nN;$j++){
	$i=$ns[$j];
	fprintf($file,"%f\t%f\t0.5\n",$xs[$i],$ys[$i]);
}
fclose($file);
/* generate data.pos*/
$content="7
y
cell.in
y
1 1 1
1
3
C N
0
CN.xyz
structure
map.in
";
		$file=fopen("in.disp","w");
	fprintf($file,$content);
	  $home=dirname(__FILE__);	
	shell_exec("$home/../../latgen <in.disp");
}
?>