  <?php
	  //$home=dirname(__FILE__);
	//require_once("$home/../../discript.php");
	function preFile(){
	global $usinglat,$xlen,$bond,$latx,$laty,$latz,$cell;
$bond=1.42;
		$dx=sqrt(3)*$bond;
		$dy=3*$bond;
	if(!$usinglat){
$latx=floor($xlen/$dx);
$laty=floor($ylen/$dy);
$latz=1;
	}
	$content="3
$dx
10
5
5
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
	}
?>