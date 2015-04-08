<?php
$srcHome="/home/xggong/home1/zhouy/tcscripts/src";
$projHome=dirname(__FILE__);
$projName=basename($projHome);
if($stage==1){
	$units="lj";
	$species="Ar";
	$method="greenKubo";
	for($lx=200;$lx<300;$lx*=2){
	submit("\$xlen=$lx * Unit('metal','l');",array("xlen"=>$lx* Unit('metal','l')));
	}
}
if(!$initsub){
shell_exec("cp $srcHome/submit.php $projHome;");
require_once("submit.php");
shell_exec("rm $projHome/submit.php");
}
#shell_exec("rm $projHome/submit.php;");
?>