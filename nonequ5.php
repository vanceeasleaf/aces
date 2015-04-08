<?php
/* 通过lammps的dump文件确定原子的近邻关系，找出不等价的原子个数*/

/* 读入坐标文件*/
$fp=fopen("../range","r");
for($i=0;$i<3;$i++)fgets($fp);
list($natom)=fscanf($fp,"%d");
$fdebug=fopen("debug.log","w");
fgets($fp);
list($xlo,$xhi)=fscanf($fp,"%f %f");
list($ylo,$yhi)=fscanf($fp,"%f %f");
list($zlo,$zhi)=fscanf($fp,"%f %f");
$lx=$xhi-$xlo;$ly=$yhi-$ylo;$lz=$zhi-$zlo;
fgets($fp);
for($i=0;$i<$natom;$i++){
	list($id,$type[$i],$x[$i],$y[$i])=fscanf($fp,"%d%s%f%f");
	$x[$i]*=$lx;
	$y[$i]*=$ly;
//	fprintf($fdebug,"$type[$i]\t$x[$i]\t$y[$i]\n");
}

/* 构建近邻列表*/
$neigh=array();
for($i=0;$i<$natom;$i++)$neigh[$i]=array();
for($i=0;$i<$natom;$i++){
	for($j=0;$j<$i;$j++){
		if(within($x[$i]-$x[$j],$y[$i]-$y[$j],1.7)){
			$neigh[$i][]=$j;
			$neigh[$j][]=$i;
		}
	}
}

/* 对所有邻居按幅角排序*/
for($i=0;$i<$natom;$i++)$neigh[$i]=sort_nei($i,$neigh[$i]);

/* 输出调试*/
$f=fopen("neighborlist.txt","w");
for($i=0;$i<$natom;$i++){
	fprintf($f,"%d",$i);
	$neii=$neigh[$i];
	for($k=0;$k<count($neii);$k++)fprintf($f,"\t%f",$neii[$k]);
	fprintf($f,"\n");
}
//exit();
/* 不等价的cluster个数*/
$nonequ=array();
for($i=0;$i<16;$i++)$nonequ[$i]=0;
for($i=0;$i<$natom;$i++){
	$code=0;
	for($j=0;$j<count($neigh[$i]);$j++){
		$k=$neigh[$i][$j];
		$dex=($type[$k]==1)?0:1;
		$code+=(1<<$j)*$dex;
	}
	$dex=($type[$i]==1)?0:1;
	$code+=(1<<3)*$dex;
	$nonequ[$code]++;
}
echo max($nonequ)/$natom;
function within($x,$y,$r){
	global $lx,$ly;
	$x=abs($x);
	$y=abs($y);
	$a=$x*$x+$y*$y<$r*$r;
	$b=($x-$lx)*($x-$lx)+$y*$y<$r*$r;
	$c=$x*$x+($y-$ly)*($y-$ly)<$r*$r;
	$d=($x-$lx)*($x-$lx)+($y-$ly)*($y-$ly)<$r*$r;
	return $a||$b||$c||$d;
}

function sort_nei($i,$neii){
	global $x,$y;
	$n=count($neii);
	
	/* 四个象限里是否有原子*/
	$occu=array(0,0,0,0);
	
	/* 给出各邻居的幅角*/
	$ang=array();
	for($k=0;$k<$n;$k++){
		$j=$neii[$k];
		$ang[$k]=atan2($y[$j]-$y[$i],$x[$j]-$x[$i])+pi();
		$idx=floor($ang[$k]/(pi()/2));
		
		/* 原子j占据了第N象限那个位置*/
		$occu[$idx]=$j;
	}
	//return $ang;
	if($occu[2]==0)return array($occu[0],$occu[1],$occu[3]);
	else return array($occu[2],$occu[3],$occu[1]);
	/* 排序
	$angs=array();
	for($k=0;$k<$n;$k++){
		$angs[$k]=array("id"=>neigh[$i][$k],"ang"=>$ang[$k]);
	}
	$angs=quicksort($angs,"ang");
	$result=array();
	for($k=0;$k<$n;$k++){
		$result[$k]=$angs[$k]['id'];
	}
	return $result;
	*/
}

/* 快排*/
function quicksort($str,$field){ 
	if(count($str)<=1) return $str;//如果个数不大于一，直接返回 
	$key=$str[0];//取一个值，稍后用来比较； 
	$left_arr=array(); 
	$right_arr=array(); 
	for($i=1;$i<count($str);$i++){//比$key大的放在右边，小的放在左边； 
		if($str[$i][$field]<=$key[$field]) 
			$left_arr[]=$str[$i]; 
			else 
			$right_arr[]=$str[$i]; 
	} 
	$left_arr=quicksort($left_arr);//进行递归； 
	$right_arr=quicksort($right_arr); 
	return array_merge($left_arr,array($key),$right_arr);//将左中右的值合并成一个数组； 
}
?>
