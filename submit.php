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
}else if($argv[1]=="q"){query();
}else if($argv[1]=="clean"){clean();
}else if($argv[1]=="stop"){stop();
}else die("Unkown command!\n");

function pwrite($fp,$s){   
	printf("$s");
	fprintf($fp,"$s");
}

function tabjoin(){
	  $argc = func_num_args();    #获取参数个数
        $argv = func_get_args();    #获取参数值
        $s="";
        for($i=0;$i<$argc;$i++){
        	$s.=$argv[$i];
        	if($i<$argc-1)$s.="\t";
        }
	return $s;
}

function getRatio($path){
	$fp=fopen($path,"r");
	fgets($fp);
	list($natom)=fscanf($fp,"%d");
	list($ntype)=fscanf($fp,"%d");
	$n=0;
	while($label!="Atoms"&&$n<20){
		list($label)=fscanf($fp,"%s");
		$n++;
	}
	fgets($fp);
	for($i=0;$i<$ntype;$i++)$a[$i]=0;
	while(list($null,$type)=fscanf($fp,"%d%d")){
		$a[$type-1]++;
	}
	return $a[1]/$natom;
}

function getRdf($path,$ratio){
	$fp=fopen($path,"r");
	for($i=0;$i<3;$i++)fgets($fp);
	list($null,$n)=fscanf($fp,"%d%d");
	$s=0;
	for($i=0;$i<$n;$i++){
		list($id,$r,$g11,$g12,$g21,$g22)=fscanf($fp,"%d%f%f%f%f%f");
		$g=-($ratio*$ratio*$g11+$ratio*(1-$ratio)*$g12*2+(1-$ratio)*(1-$ratio)*$g22);
		$p=$r*$r*($g-1);
		$s+=$p;
	}
	return $s;
}

function query(){
	global $projHome,$universe,$srcHome,$php;
	$result=fopen("$projHome/result.txt","w");
	
	  /**
	  * qloops.txt的处理分为以下几种情况:
	  * 正常情况：通过php sub.php投任务并生成qloops.txt以后
	  * 工程被复制以后
	  * 工程被移动以后
	  * 该文件不存在时
	  */
        $obj=getObjs("$projHome/qloops.txt");

	  /* work的固有属性*/
         echo "id\tpercent\tstatus\tqueue\tprocs";
         fprintf($result,"id");
         
         /* work覆盖过的参数*/
         $paras=getParas($obj);
         foreach ($paras as $para){
         	 pwrite($result,"\t".$para);
         }
         
         /* work的计算结果*/
         pwrite($result,"\tkappa\ttotalE\tNatom\tE/N\tdisorder\trd\tdisorderC\tratio\trdfs\n");
         
         checkUniverse($projHome,$universe);

        /* 遍历projet中的所有work*/
        foreach($obj as $ob){
        	$id=$ob["id"];
        	$pid=$ob["pid"];
        	$runTime=$ob["runTime"];
        	if(!$runTime)$runTime=10000000;
     	/* work的固有属性*/
        	list($percent,$status,$queue,$nodes,$procs)=getQueryInfo("$projHome/$id",$pid,$runTime,$ob);
        	echo tabjoin($id,$percent,$status,$queue,$nodes."x".$procs);
        	fprintf($result,"$id\t");
        	
        	/* work覆盖过的参数*/
        	for($j=0;$j<count($paras);$j++){
        		$key=$paras[$j];
	         	 printf("\t%s",$ob[$key]===""?"def":$ob[$key]);
	         	 if($percent+0>.5){
	         	    	 fprintf($result,"%s\t",$ob[$key]==""?"def":$ob[$key]);
	         	 }
         	 }
         	 
         	 /* work的计算结果*/
         	 if($percent+0>0){
         	          	  
 	          	  /* 准备参数列表并调用后处理*/
 	          	  $dir="$projHome/$id";
 	          	  $dir=preg_replace("/\//","\\\/",$dir);
 	          	  $sed="sed 's/projHome=.\+/projHome=\"".$dir."\";/g ' qloop.php>qloop.php1";
 	          	  if(file_exists("post.php"))$postfile= "../post.php";
 	          	  else $postfile="";
 	          	  passthru("cd $projHome/$id;$sed;cat qloop.php1 $postfile>qloop.php2;$php $srcHome/profile.php \"$projHome/$id/qloop.php2\" \"$projHome/$id/species.php\";");
 	          
         	          /* 取出后处理结果，热导率*/
         	          $kappaline=shell_exec("cd $projHome/$id;tail -1 result.txt;");
         	          $kappa=trim(substr($kappaline,strpos($kappaline,'=')+1));
         	          pwrite($result,"$kappa");
         	          
         	          /* 总能量*/
         	          $totalEline=shell_exec("cd $projHome/$id/minimize;tail -22 log.out| head -1;");
         	          list($null,$totalE)=sscanf($totalEline,"%d%f");
         	          pwrite($result,"\t$totalE");
         	        
         	          /* 原子数和平均能量*/
         	          $Natomline=shell_exec("cd $projHome/$id/minimize;head -5 log.out|tail -1 ;");
         	          list($Natom)=sscanf($Natomline,"%d");
         	          if($Natom){
         	        //pwrite($result,"\t$Natom");
         	          $epn=$totalE/$Natom;        	          
         	        //pwrite($result,"\t$epn");
         	          }
         	          
         	          /* 无序度*/
         	           $disorderLine=shell_exec("cd $projHome/$id/minimize;mkdir disorder 2>err;cd disorder;cp $srcHome/in.disorder .;$APP_PATH<in.disorder 2>err 1>log;tail -1 disorder.txt  2>err;");
         	          list($null,$disorder,$rd)=sscanf($disorderLine,"%d%f%f");
         	          pwrite($result,"\t$disorder\t$rd");
         	      //    $disorderLine=shell_exec("cd $projHome/$id/minimize;mkdir disorderdist 2>err;cd disorderdist;cp $srcHome/indist.disorder .;$APP_PATH<indist.disorder 2>err 1>log;tail -1 disorder.txt  2>err;");
         	     //   list($null,$disorder,$rd)=sscanf($disorderLine,"%d%f%f");
         	          //pwrite($result,"\t$disorder\t$rd");
         	         // $disorderLine=shell_exec("cd $projHome/$id/minimize;mkdir disorderC 2>err;cd disorderC;cp $srcHome/in.disorderC .;$APP_PATH<in.disorderC 2>err 1>log;tail -1 disorder.txt  2>err;");
         	          //list($null,$disorderC)=sscanf($disorderLine,"%d%f");
         	          //pwrite($result,"\t$disorderC");
         	          
         	          /* 掺杂比例*/
         	         // $ratio=getRatio("$projHome/$id/minimize/structure");
         	          pwrite($result,"\t$ratio");
         	          
         	          //$rdfs=getRdf("$projHome/$id/minimize/disorder/rdf.txt",$ratio);
         	          //pwrite($result,"\t$rdfs");
         	          /*
         	          $nonequ=shell_exec("cd $projHome/$id/minimize;mkdir nonequ 2>err;cd nonequ;$php $srcHome/nonequ.php;");
         	          pwrite($result,"\t$nonequ");
         	          $nonequ3=shell_exec("cd $projHome/$id/minimize/nonequ ;$php $srcHome/nonequ3.php;");
         	          pwrite($result,"\t$nonequ3");
         	           $nonequ4=shell_exec("cd $projHome/$id/minimize/nonequ ;$php $srcHome/nonequ4.php;");
         	          pwrite($result,"\t$nonequ4");*/
         	          $species=$ob["species"];
         	          if(in_array($species ,array("CN-small"))){
         	          $nonequ5=shell_exec("cd $projHome/$id/minimize;mkdir nonequ 2>err;cd nonequ;python $srcHome/inequality.py;");
         	          pwrite($result,"\t$nonequ5");
         	          }
         	          }
         	       pwrite($result,"\n");
       }
}

/**
 * 为多任务并轨准备文件
 * @author zhouy
 */
function checkUniverse($projHome,$universe){
	if(!$universe)return;
	$ff=fopen("$projHome/pbs/info","r");
	$n=0;
	while(list($uid[$n],$upid[$n],$ulog[$n],$uscreen[$n])=fscanf($ff,"%d%s%s%s")){
		$n++;
	}
	for($i=0;$i<count($obj);$i++){
		$id=$ob["id"];
		$ob["pid"]=intval($upid[$i]);
		if(!is_file("$projHome/pbs/$ulog[$i]"))continue;
		if(!is_file("$projHome/pbs/$uscreen[$i]"))continue;
		shell_exec("cp $projHome/pbs/$ulog[$i] $projHome/$id/log.lammps");
		shell_exec("cp $projHome/pbs/$uscreen[$i] $projHome/$id/log.out");
		shell_exec("cp $projHome/pbs/minimize/$ulog[$i] $projHome/$id/minimize/log.lammps");
		shell_exec("cp $projHome/pbs/minimize/$uscreen[$i] $projHome/$id/minimize/log.out");
	}
}

/**
 * 获取任务的固有属性
 * @author zhouy
 */
function getQueryInfo($workPath,$pid,$runTime,$ob){
	$lastline=@shell_exec("tail -1 $workPath/log.out");
	$qstat=@shell_exec("qstat $pid 2>&1|tail -1 ");
	list($step)=sscanf($lastline,"%d");
	$percent=sprintf("%.1f%%",$step/$runTime*100);
	if(strpos($qstat,"Unknown Job Id")){/*该任务已从任务队列中去除*/
		$time="complete";
		if(strpos($lastline,"builds")){
			$status="C";
			$percent="100%";
		}else{/*异常退出*/
			$status="E";
		}
		$queue=$ob["queue"];
		$nodes=$ob["nodes"];
		$procs=$ob["procs"];
	}else{/*正在运行或等待R&Q&C*/
		list($null,$null,$null,$time,$status,$queue)=sscanf($qstat,"%s%s%s%s%s%s");
		$info=shell_exec("qstat -f $pid 2>&1|grep nodes");
		list($null,$null,$info)=sscanf($info,"%s%s%s");
		$nnn=split(":ppn=",$info);
		$nodes=$nnn[0];
		$procs=$nnn[1];
	}
	return array($percent,$status,$queue,$nodes,$procs);
}

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
		/*
//复制文件夹后容易kill掉原来工程的程序
		$qloop=fopen("qloops.txt","r");
$n=0;			
while($json_string=fgets($qloop)){
        $obj[$n++]=json_decode($json_string,true);
        }
        for($i=0;$i<count($obj);$i++){
        	$ob=$obj[$i];
        	$pid=$ob["pid"];
	shell_exec("qdel $pid 2>log");
        }
        */
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

/**
 * 获取有效参数，排除那些已经考虑过的
 * @author zhouy
 */
function getParas($obj){
	$paras=array();
	foreach($obj as $pa){
		while($key=key($pa)){
		 if($key=="id"||
		 $key=="pid"||
		 $key=="time"||
		 $key=="cmd"||
		 $key=="project"||
		 $key=="nodes"||
		 $key=="procs"||
		 $key=="species"||
		 $key=="units"||
		 $key=="method"){next($pa);continue;}
		 $flag=0;
		 foreach($paras as $para){
		 	 if($key==$para){$flag=1;break;}
		 }
		 if($flag==0)array_push($paras,$key);
		next($pa);
		}
     }
     return $paras;
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
