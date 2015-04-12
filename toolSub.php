<?php
$idx=0;

/**
 * 获取参数并生成qloops.php 将参数保存，生成lammps.pbs等队列信息，并投任务
 * @author zhouy
 */
function submit(){
	global $projHome;
	global $stage;
	if($stage!=1)return;
	global $projName;
	$argc= func_num_args();    #获取参数个数
       $argv = func_get_args();    #获取参数值
       if($argc>0)$cmd=$argv[0];
       else $cmd="";
	global $idx;
	global $loops;
	global $nodes;
	global $procs;
	global $species;
	global $method;
	global $units;
	global $universe;
	global $queue;
	global $runTime,$uqueue,$single,$unodes,$uprocs;
	global $srcHome;
	if(!$queue)$queue="q1.1";
	if(!$nodes)$nodes=1;
	if(!$procs)$procs=12;
	$jj=json_encode($argv);
	passthru("python $srcHome/toolsub.py $projHome $projName '$cmd' $idx $nodes $procs $species $method $units '$universe' $queue $runTime '$uqueue' '$single' '$unodes' '$uprocs' '$jj'");
	$idx++;
}

/**
 * 统一submit的两个参数，只提供它的第二个参数，自动生成第一个参数
 * @author zhouy
 */
function submitq($pa=array()){
	$cmd="";
	foreach($pa as $key=>$val){
		$cmd.="\$".$key."=$val;";
	}
	submit($cmd,$pa);
}

/**
 * 将队列的多个核拆分为几组并运行几个work
 * @author zhouy
 */
function uexec(){
	global $universe;
	if(!$universe)return;
	global $projHome;
	$home=dirname(__FILE__);
	shell_exec("cd $projHome/pbs;ls *.pbs>tmp;");
	$fp=fopen("$projHome/pbs/tmp","r");
	$n=0;
	while(list($st[$n],$ed[$n])=fscanf($fp,"lammps%d-%d.pbs")){
		$n++;
	}
	
	/* sort*/
	$len=$n;
	for($i=0;$i<$len;$i++){
		$idx[$i]=$i;
	}
	for($i=0;$i<$len;$i++)
		for($j=0;$j<$i;$j++){
			if($st[$i]<$st[$j]){
				swap($st[$i],$st[$j]);
				swap($ed[$i],$ed[$j]);
				swap($idx[$i],$idx[$j]);
			}
		}
	$fi=fopen("$projHome/pbs/info","w");
	for($i=0;$i<$len;$i++){
		$p=$i;
		$st1=$st[$p];
		$ed1=$ed[$p];
		echo "lammps${st1}-${ed1}.pbs\n";
		$lb="${st1}-${ed1}";
		$pid=shell_exec("cd $projHome/pbs;qsub lammps${st1}-${ed1}.pbs;");
		echo $pid;
		$pid=trim($pid);
		for($j=$st1;$j<=$ed1;$j++){
			$b=$j-$st1;
			fprintf($fi,"$j\t$pid\tlog.$lb.$b\tscreen.$lb.$b\n");
		}
	}	
}	

/* exchange*/
function swap(&$a,&$b){
	$tmp=$b;$b=$a;$a=$tmp;
}

?>
