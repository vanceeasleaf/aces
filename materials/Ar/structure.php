  <?php
	  $home=dirname(__FILE__);
	require_once("$home/../../discript.php");
	$file=fopen("in.disp","w");
	if(!$usinglat){
$latx=floor($xlen/(2*$bond));
$laty=floor($ylen/(2*$bond));
$latz=floor($thick/(2*$bond));
	}
	echo "
lattice            fcc    $latticeConstant
region            box  block 0  $latx  0 $laty  0 $latz  units lattice
create_box     1  box
create_atoms     1  box
";
?>