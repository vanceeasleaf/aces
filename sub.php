<?php
$srcHome="/home/xggong/home1/zhouy/tcscripts/src";
$projHome=dirname(__FILE__);
$projName=basename($projHome);

	$units="metal";
	$species="SiGe";
	$method="greenkubo";
	$nodes=1;
	$procs=4;$queue="q1.4";
	
	$runTime=10000000;
	if($stage==1){
	for($i=0;$i<20;$i++){
	submit("\$seed=13513;\$computeTc=1;\$cell=\"SiGesuper/$i\";\$usinglat=1;\$hdeta=\$deta;\$latx=1;\$laty=1;\$latz=1;");
	}
}
shell_exec("cp $projHome/sub.php $srcHome;");
require_once("$srcHome/submit.php");
?>
