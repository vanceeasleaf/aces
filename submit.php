<?php
$initsub=1;
require_once("sub.php");
require_once("$srcHome/config.php");
require_once("$srcHome/funcs.php");
$home=dirname(__FILE__);
require_once("$srcHome/exec.php");

if($argc==1){
	passthru("cd $projHome;python $srcHome/query.py clean $projHome $srcHome '$universe' $php $projName '$single'");
	require_once("$srcHome/toolSub.php");
	$stage=1;
	require("sub.php");
	exit();
}else if($argv[1]=="q"){//query();
}else if($argv[1]=="clean"){//clean();
}else if($argv[1]=="stop"){//stop();
}else die("Unkown command!\n");
$aa=$argv[1];
passthru("cd $projHome;python $srcHome/query.py $aa $projHome $srcHome '$universe' $php $projName '$single'");
?>
