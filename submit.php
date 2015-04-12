<?php
require_once("sub.php");
$home=dirname(__FILE__);
if($argc==1){
	//clean
	passthru("cd $projHome;python $srcHome/query.py clean $projHome $srcHome '$universe' $projName '$single'");
	require_once("$srcHome/toolSub.php");
	$stage=1;
	require("sub.php");
	exit();
}else if($argv[1]=="q"){//query();
}else if($argv[1]=="clean"){//clean();
}else if($argv[1]=="stop"){//stop();
}else die("Unkown command!\n");
$aa=$argv[1];
passthru("cd $projHome;python $srcHome/query.py $aa $projHome $srcHome '$universe' $projName '$single'");
?>
