<?php
$srcHome="/home/xggong/home1/zhouy/tcscripts/src";
$projHome=dirname(__FILE__);
$projName=basename($projHome);
if($stage==1){
	$units="metal";
	$species="CN-small";
	$method="nvt";
	$nodes=1;
$procs=4;$queue="q1.4";
$runTime=10000000;
for($i=0;$i<31;$i++)	
submit("\$cell=\"C51N/$i\";\$thick=1.44;\$langevin=0;\$hdeta=8*\$deta;\$usinglat=1;\$timestep=.5e-3;\n\$latx=11;\$laty=1;\$latz=1;");
}
shell_exec("cp $projHome/sub.php $srcHome;");
require_once("$srcHome/submit.php");
?>
