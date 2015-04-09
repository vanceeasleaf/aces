<?php
$initsub=1;
require_once("sub.php");
require_once("$srcHome/config.php");
require_once("$srcHome/funcs.php");
$home=dirname(__FILE__);
require_once("$srcHome/exec.php");

if($argc==1){
	clean();
	require_once("$srcHome/toolSub.php");
	$stage=1;
	require("sub.php");
}else if($argv[1]=="q"){//query();
passthru("cd $projHome;python $srcHome/query.py $projHome $srcHome '$universe' $php");
}else if($argv[1]=="clean"){clean();
}else if($argv[1]=="stop"){stop();
}else die("Unkown command!\n");


function stop(){
	printf("Comfirm to stop all the simulation in this project?[y/n]");
	$stdin = fopen('php://stdin', 'r');
	list($s)=fscanf($stdin,"%s"); 
	if($s!="y")exit("exit with no change.");
	comfirmStop();
}

function comfirmStop(){
	global $projName;
	global $projHome;
	global $single;
	
	/* 容易kill掉同名工程程序*/
	if($single){
		$obj=getObjs("$projHome/qloops.txt");
	      for($i=0;$i<count($obj);$i++){
         	 	$pa=$obj[$i];
         	 	$pid=$pa["pid"];
         	 	echo "kill:$pid\n";
         	 	exec::kill($pid);
         	 }
	return;
	}
	$tarname="zy_$projName"."_";
	if(strlen($tarname)<10){
		shell_exec("qstat|grep $tarname>tmp");
		$file=fopen("$projHome/tmp","r");
		while(list($pid)=fscanf($file,"%s")){
			$pid=intval($pid);
			echo "qdel:$pid\n";
			shell_exec("qdel $pid 2>log");
		}
	}else{
		shell_exec("qstat|grep xggong >tmp");
		$file=fopen("$projHome/tmp","r");
		while(list($pid)=fscanf($file,"%s")){
			$pid=intval($pid);
			$jobnameString=shell_exec("qstat -f $pid |grep Job_Name");
			list($null,$null,$jobname)=sscanf($jobnameString,"%s%s%s");
			if(strstr($jobname,$tarname)){
				echo "qdel:$pid\n";
				shell_exec("qdel $pid 2>log");
			}
		}
	}

function getObjs($fileName){
	
	/*执行程序之前会清理上次的，读取qloop.txt，如果没有上次的文件就不清理。*/
	if(!is_file($fileName))return;
		$qloop=fopen($fileName,"r");
		$n=0;
		while($json_string=fgets($qloop)){
        	$obj[$n++]=json_decode($json_string,true);
        }
        return $obj;
}

function clean(){
	global $projHome;
	printf("Comfirm to clean all the files in this project?[y/n]");
	$stdin = fopen('php://stdin', 'r');
	list($s)=fscanf($stdin,"%s"); 
	if($s!="y")exit();
	comfirmStop();
	
	/*删除原始代码以外的文件*/	
	$files=shell_exec("cd $projHome;ls ");
	$files=explode("\n",$files);
	foreach ($files as $ls){
		if($ls=="sub.php"||$ls=="post.php"||$ls=="data"||$ls=="")continue;
		echo "deleting:$ls\n";
		shell_exec("cd $projHome;rm -r $ls");
	}
}

?>
