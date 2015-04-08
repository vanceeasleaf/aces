  <?php
	  $home=dirname(__FILE__);
	require_once("$home/../../discript.php");
$bond=1.42*Unit($srcUnit,"l");
		$dx=sqrt(3)*2*$bond;
		$ax=sqrt(3)*2;
		$dy=6*$bond;
		$ay=6;
	if(!$usinglat){
$latx=floor($xlen/$dx);
$laty=floor($ylen/$dy);
$latz=1;
	}
		$file=fopen("cell.in","w");
$content="1
$dx 0 0 
0 $dy 0
0 0 100
12 4
";
$content=sprintf("$content%f\t%f\t%f\n",0.000000,1.0/$ay,0.000000);
$content=sprintf("$content%f\t%f\t%f\n",0.866050/$ax,1.500043/$ay,0.000000);
$content=sprintf("$content%f\t%f\t%f\n",0.866050/$ax,    2.500071/$ay,    0.000000);
$content=sprintf("$content%f\t%f\t%f\n",1.732100/$ax,    0.000000/$ay,    0.000000);
$content=sprintf("$content%f\t%f\t%f\n",0.866050/$ax,    4.500128/$ay,    0.000000);
$content=sprintf("$content%f\t%f\t%f\n",0.866050/$ax,    5.500156/$ay,    0.000000);
$content=sprintf("$content%f\t%f\t%f\n",1.732100/$ax,    4.000113/$ay,    0.000000);
$content=sprintf("$content%f\t%f\t%f\n",2.598150/$ax,    4.500128/$ay,    0.000000);
$content=sprintf("$content%f\t%f\t%f\n",2.598150/$ax,    5.500156/$ay,    0.000000);
$content=sprintf("$content%f\t%f\t%f\n",2.598150/$ax,    1.500043/$ay,    0.000000);
$content=sprintf("$content%f\t%f\t%f\n",2.598150/$ax,    2.500071/$ay,    0.000000);
$content=sprintf("$content%f\t%f\t%f\n",0.000000/$ax,    3.000085/$ay,    0.000000);
$content=sprintf("$content%f\t%f\t%f\n",0.000000/$ax,    0.000000/$ay,    0.000000);
$content=sprintf("$content%f\t%f\t%f\n",1.732100/$ax,    1.000028/$ay,    0.000000);
$content=sprintf("$content%f\t%f\t%f\n",0.000000/$ax,    4.000113/$ay,    0.000000);
$content=sprintf("$content%f\t%f\t%f\n",1.732100/$ax,    3.000085/$ay,    0.000000);
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