  <?php
	  //$home=dirname(__FILE__);
	//require_once("$home/../../discript.php");
function rotate($phi,$x,$y){
	$x1=cos($phi)*$x-sin($phi)*$y;
	$y1=sin($phi)*$x+cos($phi)*$y;
	return array($x1,$y1);
}
function preFile(){
	global $usinglat,$xlen,$bond,$latx,$laty,$latz,$cell;
$pos1=array(
6.928400,13.000369,0.000000
	,7.794450,16.500469,0.000000
	);
$phi=3.14159265359/2-atan(($pos1[4]-$pos1[1])/($pos1[3]-$pos1[0]));
$bond=sqrt(($pos1[4]-$pos1[1])*($pos1[4]-$pos1[1])+($pos1[3]-$pos1[0])*($pos1[3]-$pos1[0]))*1.42;
		$dx=sqrt(3)*$bond*2;
		$dy=3*$bond*2;
	if(!$usinglat){
$latx=floor($xlen/$dx);
$laty=floor($ylen/$dy);
$latz=1;
	}
$home=dirname(__FILE__);
shell_exec("cp $home/pos/$cell cell.in");
$content="7
y
cell.in
y
$latx $laty $latz 
1
3
C C
0
CD.xyz
structure
map.in
";
		$file=fopen("in.disp","w");
	fprintf($file,$content);
	  $home=dirname(__FILE__);	
	shell_exec("$home/../../latgen <in.disp");
}
?>