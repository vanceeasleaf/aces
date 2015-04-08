  <?php
	  //$home=dirname(__FILE__);
	//require_once("$home/../../discript.php");
function rotate($phi,$x,$y){
	$x1=cos($phi)*$x-sin($phi)*$y;
	$y1=sin($phi)*$x+cos($phi)*$y;
	return array($x1,$y1);
}
function preFile(){
	global $usinglat,$xlen,$bond,$latx,$laty,$latz,$cell,$ratio;

	if(!$usinglat){
$latx=floor($xlen/$bond);
$laty=floor($ylen/$bond);
$latz=floor($thick/$bond);
	}
$home=dirname(__FILE__);

if($cell){
	shell_exec("cp $home/pos/$cell cell.in");
	$content=
"7
y
cell.in
y
$latx $laty $latz 
1
3
Si Ge
0
SiGe.xyz
structure
map.in
";
}else{
$content=
"1
$bond
6
1
y
$latx $laty $latz 
1\n"//x go fast
."1\n"//substituional solid solution
."5\n"//all
."1\n"//the type to be replaced
."$ratio\n"
."2\n"//new type
."3
Si Ge
0
SiGe.xyz
structure
map.in
";
}
		$file=fopen("in.disp","w");
	fprintf($file,$content);
	  $home=dirname(__FILE__);	
	shell_exec("$home/../../latgen <in.disp");
}
?>