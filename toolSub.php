<?php
require_once("genPbs.php");
require_once("genPbss.php");
require_once("units.php");
require_once("$srcHome/exec.php");
if(!$projHome)exit("please set your project path!");

/* 设置时区*/
date_default_timezone_set("PRC");

$loops=fopen("$projHome/qloops.txt","w");
$idx=0;

/**
 * 运行某个任务
 * @author zhouy
 * @input 任务名称
 * @output 任务队列号
 */
function setSubProject($index){
	global $projHome;
	global $single;
	if($single)$pid=exec::background("sh $projHome/$index/run.sh");//$pid=exec("cd $projHome/$index;sh run.sh > /dev/null & echo $!");
	else
	$pid=shell_exec("cd $projHome/$index;qsub lammps.pbs;");
	echo "submit: $pid\t$projHome/$index\n";
	#sleep(1);
	return intval($pid);
}

function write($cmd,$fileName){
	$file=fopen($fileName,"w");
	fprintf($file,$cmd."\n");
	fclose($file);
}

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
	global $runTime;
	if(!$queue)$queue="q1.1";
	if(!$nodes)$nodes=1;
	if(!$procs)$procs=12;
	$para="";
	
	/* 创建任务序列文件，并把submit传参过来的交给qloop.php*/
	makeLoopFile($cmd,$idx);
	if(!$universe)$pid=setSubProject($idx);
	
	/* project中每个work的信息需要保留给后续进程*/
	$json_obj=array(
		"id"=>$idx,
		"pid"=>$pid,
		"time"=>date('Y-m-d H:i:s'),
		"cmd"=>$cmd,
		"nodes"=>$nodes,
		"procs"=>$procs,
		"species"=>$species,
		"method"=>$method,
		"units"=>$units,
		"runTime"=>$runTime
	);
		
	/* 使用submit的参数传递的参数*/
	if($argc>1){
		$pa=$argv[1];
		while($key=key($pa)){
				$json_obj[$key]=$pa[$key];
				next($pa);
		}
	}
	
	/* 每行一个json*/
	fprintf($loops,"%s\n",json_encode($json_obj));
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

/**
 * 生成某个work的参数列表
 * @author zhouy
 */
function makeLoopFile($cmd,$idx){
	global $projHome;
	global $projName;
	global $species;
	global $units;
	global $method;
	global $queue;
	global $nodes;
	global $procs;
	global $uqueue;
	global $single;
	global $universe;
	
	/**
	 * universe ,single ,normal是几种不同的控制方式 ,应该有一个接口描述他们。project掌管高层逻辑，不想知道任务是怎么分配的；而work只管自己的输入和结果，对
	 * 任务细节不感兴趣。
	 * 那么系统的最高级逻辑为：
	 * 用户想要执行一些work，分别指定参数和类型，希望有一个指挥系统，可以统一调度这些work，也可以选择性地调度其中几个。
	 * 调度操作包括 work的前处理，work的执行，work的后处理，work的中止，目前不需要选择性调度。
	 * 可以比喻为一个多功能烤箱，用户在里面放各种饼干，烤箱给他们刷油（前处理），烘焙（执行），脱模（后处理），当然也可以随时中止。同时有一些反馈，比如亮灯。
	 * 但并不在乎刷油是一个接一个地刷还是用大刷子一行一行地刷。
	 * 故project实现了Userinterface接口。project拥有一个controller对象，用于决定并行或串行等执行方式（用户可选）。
	 * 现在考虑另一个需求：烘焙过程中查看饼干的进度。这种饼干一旦烤好了，再怎么加热都不会再变，即100%烤好。每个饼干烤好所需时间不同，问烤箱如何计算这种进度。
	 * 采用问的方式即可：烤箱问饼干，现在是20步，你的进度是多少。即增加饼干的API。但实际上烤箱怎么知道到第几步了，它只管加热，到第几步跟饼干的吸收还有关系。
	 * 因此应该直接问饼干，你的进度是多少。
	 * 最后的反馈并非所有饼干烤好闪两下灯就行了，还得由烤箱生成一份此次烘焙的饼干的报告。当然一切的信息都应该是问饼干得到的，至于哪些信息应该问，不如把控制本次
	 * 烘焙的所有更改过参数都问了吧，以及饼干烘焙好后的测量结果。结果打印出来备案。
	 * 烤箱的屏幕需要输出缩略信息。
	 * 烘焙后的测量结果有哪些要问由用户决定。为了结果的可比性，我们要求这些测量结果每个饼干都应该具有，要不然他们不能放在一起烤。
	 * 烤箱事先要问饼干有哪些结果是可以提供的。
	 */
	 
 	  /**
	   * 现在剩下一个问题，对于烤箱和饼干都需要的参数应该如何处理。具体指projHome等参数。生成各个work的文件结构的工作是由谁来做，work为什么需要projHome等工程的自有
	   * 参数？是由projct类来做？这部分工作属于前处理，但是project作为一个指挥系统当然不知道如何进行work的前处理，因此还是要委托给work自己来做。它只是一个userinterface而已，
	   * 真正进行复杂工作的不能是它。work获取projHome参数的目的是创建目录结构，unit参数的存在是因为给work提供参数的时候要指定单位。species和method对工程而言没有作用。
	   * 总的来说真正共享的参数是projhome/idx（即work需要知道自己的目录）和units。似乎工程本身并不需要知道单位。
	   * 另一方面，queue,nodes,procs这几个参数是给controller的，根据情况还可能是uqueue,unodes,uprocs，或之后可实现的adaptive自适应，work本身应该知道它的运行参数吗？
	   * 就像是烤箱知道饼干的材质不同，所以对他们施加不同的功率，可能开一个大灯均一处理，也可能每个开一个小灯处理，同时还可以分成几个小组，每个组共用一个灯。最常见的情况其实
	   * 是每个开一个小灯，但灯的功率一样。
	   * 烤箱可以先行测量饼干的材质，然后再做打算。
	   * 灯这个对象应该存在吗? 当project预处理的时候生成了pbs文件，它全权决定了这个或这组饼干的烘焙过程。
	   * 所以当前的设计可以表述为：project构成用户界面，用户实例化多个work（指定相应参数），project前处理提交给controlller，controller给饼干分组并配置灯，
	   * 饼干做好准备（work展开自己的文件）。
	   * project执行，各灯点亮。
	   * query阶段，project询问各work他们的信息并进行整合。反馈到屏幕和文件上。
	   */
	   
	   /**
	    * 设计方面的问题都解决了。
	    * 然而程序并非驻留在内存中的，query阶段必须知道submit阶段的信息。这些信息留在qloop.txt中。应该包括project重建（工程名），controller（universe类型及其参数等），
	    * heaters（它的参数，以及包括哪些work），work的参数。总之query阶段重建了所有类并设定一些后处理参数，比如all.set("upP",3);
	    */
	    
	    /** 
	     * heater在执行的时候调用pbs，pbs安排work接下来的工作。
	     */
	if($universe){
		printf("prepared: $projHome/$idx\n");
		global $unodes;
		global $uprocs;
		$cores=$procs*$nodes;
		if(!$unodes)$unodes=20;
		if(!$uprocs)$uprocs=1;
		$ucores=$unodes*$uprocs;
		$len=floor($ucores/$cores);
		if(!$uqueue)$uqueue="q3.4";
		if($idx%$len==0)
		genPbss("$projHome","zy_$projName"."_",$uqueue,$unodes,$uprocs,$idx,$cores);
	}

	shell_exec("mkdir -p $projHome/$idx;cd $projHome/$idx;mkdir -p minimize");
	if($single){
		genSh("$projHome/$idx","zy_$projName"."_$idx",$procs);
	}
	if(!$universe&&!$single)genPbs("$projHome/$idx","zy_$projName"."_$idx",$queue,$nodes,$procs);
	
	/* 通过submit传递的参数*/
	write("<?php\n$cmd;\n\$projHome=\"$projHome/$idx\";\n?>","$projHome/$idx/qloop.php");
	
	/* 通过sub.php传递的参数*/
	write("<?php\n\$species=\"$species\";\n\$units=\"$units\";\n\$method=\"$method\";\n?>","$projHome/$idx/species.php");
	
	/**
	 * 传参系统的统一讨论：
	 * 当前系统中存在这样几种传参： sub.php 的变量传参->工程系统（如$method)，submit中传参-> 执行系统(如submit("\$runTime=1000000");)，
	 * sub.php变量传参->执行系统（species.php) ，submit的第二个传参->json
	 * sub.php+submit->后续系统（query目前使用的方案是：将submit过程的传参重复一遍，即以上3种传参，当后处理需要使用application中的参数时就无能为力
	 * 如runTime本来只需要给执行系统，但是为了计算进度就需要在sub.php中指定并->json，且不能靠submit的第二个传参，因为如果用户想要使用默认值的话无法获取。
	 * submit的两套传参可以统一用 submitq,但它依然无法获取默认值。即： 工程系统希望知道执行系统的参数是不可能的。
	 * sub.php向工程系统传参的目的是$nodes等pbs相关的参数。
	 * 且submit提供的参数执行在后面，还是必须要species.php
	 */
	 
	 /**
	  * 现在我们要复用参数系统来实现shengbte,并且要增加几种计算热导率的方法，增加几个需要计算的物理量。
	  */
	  

}
?>
