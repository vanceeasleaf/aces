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
$pos2=array(
0.34109,1.18154,0
,7.16273,4.33233,0
,4.77515,3.74156,0
,5.11624,2.36309,0
,2.38758,3.15078,0
,5.45733,6.10465,0
,7.50382,2.95387,0
,3.41082,4.13541,0
,2.72866,1.77231,0
,0.68216,4.92311,0
,8.18599,5.31696,0
,4.09299,1.37846,0
,8.18599,0.19692,0
,6.82166,0.59078,0
,2.04649,14.76931,0
,3.41082,14.37546,0
,5.45733,11.22468,0
,1.36433,12.40622,0
,0.68216,10.04313,0
,0.34109,11.42161,0
,6.13949,13.58776,0
,3.06975,10.6339,0
,2.72866,12.01237,0
,7.8449,11.81546,0
,5.11624,12.60315,0
,7.50382,13.19392,0
,8.52706,9.0585,0
,6.82166,10.83083,0
,6.48057,7.08926,0
,3.75191,7.87696,0
,4.77515,8.86159,0
,6.13949,8.46774,0
,1.70542,5.90772,0
,1.36433,7.2862,0
,4.09299,6.4985,0
,2.04649,9.64928,0
,0,2.56002,0
,0,7.68004,0
,9.8914,3.54463,0
,14.66656,4.72618,0
,13.6433,3.74156,0
,14.32548,6.10465,0
,16.37197,2.95387,0
,11.9379,5.51387,0
,12.27897,4.13541,0
,11.59681,1.77231,0
,9.55032,4.92311,0
,17.05414,5.31696,0
,16.71305,6.69542,0
,12.96114,1.37846,0
,10.57357,0.7877,0
,17.05414,0.19692,0
,15.34872,1.96924,0
,15.68981,0.59078,0
,10.91464,14.76931,0
,17.39521,14.17854,0
,12.27897,14.37546,0
,14.32548,11.22468,0
,12.62006,12.997,0
,10.23248,12.40622,0
,9.8914,13.7847,0
,9.55032,10.04313,0
,15.00764,13.58776,0
,11.59681,12.01237,0
,14.66656,14.96624,0
,16.37197,13.19392,0
,15.68981,10.83083,0
,13.6433,8.86159,0
,15.00764,8.46774,0
,16.03088,9.45235,0
,10.23248,7.2862,0
,12.96114,6.4985,0
,11.25573,8.27081,0
,10.91464,9.64928,0
,8.86815,2.56002,0
,8.86815,7.68004,0
,13.30223,0,0
,13.30223,10.24006,0
,1.02325,18.90472,0
,5.7984,20.08626,0
,4.77515,19.10165,0
,5.45733,21.46474,0
,7.50382,18.31395,0
,3.06975,20.87396,0
,3.41082,19.4955,0
,2.72866,17.1324,0
,0.68216,20.28319,0
,8.18599,20.67704,0
,7.8449,22.0555,0
,4.09299,16.73855,0
,1.70542,16.14778,0
,8.18599,15.557,0
,6.48057,17.32933,0
,6.82166,15.95087,0
,2.04649,30.12939,0
,8.52706,29.53863,0
,3.41082,29.73555,0
,5.45733,26.58476,0
,3.75191,28.35709,0
,1.36433,27.76631,0
,1.02325,29.14478,0
,0.68216,25.40322,0
,6.13949,28.94785,0
,2.72866,27.37246,0
,5.7984,30.32633,0
,7.50382,28.554,0
,6.82166,26.19091,0
,4.77515,24.22167,0
,6.13949,23.82783,0
,7.16273,24.81244,0
,1.36433,22.64628,0
,4.09299,21.85859,0
,2.38758,23.6309,0
,2.04649,25.00937,0
,0,17.92011,0
,0,23.04013,0
,4.43408,15.36009,0
,4.43408,25.60015,0
,9.20924,16.54163,0
,16.03088,19.69242,0
,13.6433,19.10165,0
,13.98439,17.72318,0
,11.25573,18.51087,0
,14.32548,21.46474,0
,16.37197,18.31395,0
,12.27897,19.4955,0
,11.59681,17.1324,0
,9.55032,20.28319,0
,17.05414,20.67704,0
,12.96114,16.73855,0
,17.05414,15.557,0
,15.68981,15.95087,0
,10.91464,30.12939,0
,12.27897,29.73555,0
,14.32548,26.58476,0
,10.23248,27.76631,0
,9.55032,25.40322,0
,9.20924,26.78169,0
,15.00764,28.94785,0
,11.9379,25.99398,0
,11.59681,27.37246,0
,16.71305,27.17554,0
,13.98439,27.96324,0
,16.37197,28.554,0
,17.39521,24.41859,0
,15.68981,26.19091,0
,15.34872,22.44935,0
,12.62006,23.23705,0
,13.6433,24.22167,0
,15.00764,23.82783,0
,10.57357,21.26781,0
,10.23248,22.64628,0
,12.96114,21.85859,0
,10.91464,25.00937,0
,8.86815,17.92011,0
,8.86815,23.04013,0
,1.02325,3.54463,0
,5.7984,4.72618,0
,3.06975,5.51387,0
,7.8449,6.69542,0
,1.70542,0.7877,0
,6.48057,1.96924,0
,8.52706,14.17854,0
,3.75191,12.997,0
,1.02325,13.7847,0
,5.7984,14.96624,0
,7.16273,9.45235,0
,2.38758,8.27081,0
,4.43408,0,0
,4.43408,10.24006,0
,9.20924,1.18154,0
,16.03088,4.33233,0
,13.98439,2.36309,0
,11.25573,3.15078,0
,9.20924,11.42161,0
,11.9379,10.6339,0
,16.71305,11.81546,0
,13.98439,12.60315,0
,17.39521,9.0585,0
,15.34872,7.08926,0
,12.62006,7.87696,0
,10.57357,5.90772,0
,0.34109,16.54163,0
,7.16273,19.69242,0
,5.11624,17.72318,0
,2.38758,18.51087,0
,0.34109,26.78169,0
,3.06975,25.99398,0
,7.8449,27.17554,0
,5.11624,27.96324,0
,8.52706,24.41859,0
,6.48057,22.44935,0
,3.75191,23.23705,0
,1.70542,21.26781,0
,9.8914,18.90472,0
,14.66656,20.08626,0
,11.9379,20.87396,0
,16.71305,22.0555,0
,10.57357,16.14778,0
,15.34872,17.32933,0
,17.39521,29.53863,0
,12.62006,28.35709,0
,9.8914,29.14478,0
,14.66656,30.32633,0
,16.03088,24.81244,0
,11.25573,23.6309,0
,13.30223,15.36009,0
,13.30223,25.60015,0
);
for($i=0;$i<208;$i++){
	$rpos[$i]=rotate(0,$pos2[$i*3],$pos2[$i*3+1]);
}
$minx=100000;
$miny=100000;
for($i=0;$i<208;$i++){
	$minx=min($rpos[$i][0],$minx);
	$miny=min($rpos[$i][1],$miny);
}
for($i=0;$i<208;$i++){
	$rpos[$i][0]-=$minx;
	$rpos[$i][1]-=$miny;
}
		$file=fopen("cell.in","w");
$content="1
$dx 0 0 
0 $dy 0
0 0 100
156 52
";
for($i=0;$i<208;$i++){
$content=sprintf("$content%f\t%f\t%f\n",$rpos[$i][0]/$dx,$rpos[$i][1]/$dy,0.000000);
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