  <?php
	  $home=dirname(__FILE__);
	require_once("$home/../../discript.php");
function rotate($phi,$x,$y){
	$x1=cos($phi)*$x-sin($phi)*$y;
	$y1=sin($phi)*$x+cos($phi)*$y;
	return array($x1,$y1);
}
$pos1=array(
6.928400,13.000369,0.000000
	,7.794450,16.500469,0.000000
	,9.526550,10.500299,0.000000
	,11.258650,17.500498,0.000000
	);
$phi=3.14159265359/2-atan(($pos1[4]-$pos1[1])/($pos1[3]-$pos1[0]));
$bond=sqrt(($pos1[4]-$pos1[1])*($pos1[4]-$pos1[1])+($pos1[3]-$pos1[0])*($pos1[3]-$pos1[0]))*1.42;
		$dx=sqrt(3)*$bond;
		$dy=3*$bond;
	if(!$usinglat){
$latx=floor($xlen/$dx);
$laty=floor($ylen/$dy);
$latz=1;
	}
$pos2=array(
7.794450,13.500383,0.000000
,6.928400,12.000340,0.000000
,12.124700,13.000369,0.000000
,11.258650,13.500383,0.000000
,10.392600,13.000369,0.000000
,10.392600,12.000340,0.000000
,8.660500,13.000369,0.000000
,11.258650,14.500412,0.000000
,12.124700,12.000340,0.000000
,9.526550,14.500412,0.000000
,9.526550,13.500383,0.000000
,8.660500,12.000340,0.000000
,7.794450,14.500412,0.000000
,12.990750,13.500383,0.000000
,12.990750,14.500412,0.000000
,9.526550,11.500326,0.000000
,7.794450,11.500326,0.000000
,12.124700,10.000284,0.000000
,11.258650,11.500326,0.000000
,11.258650,10.500299,0.000000
,10.392600,21.000597,0.000000
,14.722850,19.500553,0.000000
,11.258650,20.500582,0.000000
,12.124700,18.000511,0.000000
,11.258650,19.500553,0.000000
,9.526550,19.500553,0.000000
,9.526550,20.500582,0.000000
,8.660500,18.000511,0.000000
,8.660500,19.000540,0.000000
,12.990750,19.500553,0.000000
,10.392600,18.000511,0.000000
,10.392600,19.000540,0.000000
,12.990750,20.500582,0.000000
,13.856800,18.000511,0.000000
,12.124700,19.000540,0.000000
,13.856800,19.000540,0.000000
,13.856800,16.000454,0.000000
,12.990750,17.500498,0.000000
,12.124700,15.000426,0.000000
,10.392600,16.000454,0.000000
,11.258650,16.500469,0.000000
,12.124700,16.000454,0.000000
,12.990750,16.500469,0.000000
,8.660500,15.000426,0.000000
,8.660500,16.000454,0.000000
,10.392600,15.000426,0.000000
,9.526550,16.500469,0.000000
,9.526550,17.500498,0.000000
	,6.928400,13.000369,0.000000
	,7.794450,16.500469,0.000000

	,11.258650,17.500498,0.000000
			,9.526550,10.500299,0.000000
);
for($i=0;$i<52;$i++){
	$rpos[$i]=rotate($phi,$pos2[$i*3],$pos2[$i*3+1]);
}
$minx=100000;
$miny=100000;
for($i=0;$i<52;$i++){
	$minx=min($rpos[$i][0],$minx);
	$miny=min($rpos[$i][1],$miny);
}
for($i=0;$i<52;$i++){
	$rpos[$i][0]-=$minx;
	$rpos[$i][1]-=$miny;
}
		$file=fopen("cell.in","w");
$content="1
$dx 0 0 
0 $dy 0
0 0 100
51 1
";
for($i=0;$i<52;$i++){
$content=sprintf("$content%f\t%f\t%f\n",$rpos[$i][0]*1.42/$dx,$rpos[$i][1]*1.42/$dy,0.000000);
}
	fprintf($file,$content);
	$content="7
y
cell.in
y
$latx $laty $latz 
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
?>