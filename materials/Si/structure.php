  <?php
	  $home=dirname(__FILE__);
	require_once("$home/../../discript.php");
	$file=fopen("in.disp","w");
	if(!$usinglat){
$latx=floor($xlen/$bond);
$laty=floor($ylen/$bond);
$latz=floor($thick/$bond);
	}
	fprintf($file,
"1
$bond
6
1
y
$latx $laty $latz 
1
3
Si
0
Si.xyz
structure
map.in
");
	  $home=dirname(__FILE__);	
	shell_exec("$home/../../latgen <in.disp");
?>