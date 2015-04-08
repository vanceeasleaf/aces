<?php
$srcHome="/home/xggong/home1/zhouy/tcscripts/src";
$projHome=dirname(__FILE__);
$projName=basename($projHome);

	$units="metal";
	$species="C2N-h2D1";
	$method="greenkubo";
	$nodes=1;
$procs=4;$queue="q1.1";
$runTime=10000000;
if($stage==1){
		$paras=array(
		array(0),array(1),array(2),array(3),array(4)
		,array(5),array(6),array(7),array(8)
	);
	foreach( $paras as $arr ){
	list($sideLen)=$arr;
submitq(array("latx"=>2,"laty"=>2,"sideLen"=>$sideLen,"dumpxyz"=>1,"timestep"=>.5e-3));
	}
}
shell_exec("cp $projHome/sub.php $srcHome;");
require_once("$srcHome/submit.php");
?>
