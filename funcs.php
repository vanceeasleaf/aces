<?php
/** 
 * 一些工具函数
 * @author zhouy
 */

/* 产生0-1之间的随机数，包括0不包括1*/
function random01()
{
	return mt_rand() / mt_getrandmax();
}

/* 删除当前目录下的一个或多个文件或文件夹*/
function rm($file){
	$paramNum = func_num_args();    
	$params = func_get_args();    
	for($i=0;$i<$paramNum;$i++){
			$path=trim(shell_exec("pwd"));
		echo "deleting:$path/$params[$i]\n";
		shell_exec("if [ -e ".$params[$i]." ]; then rm -r ".$params[$i]."; fi");
	}
}

/* 获得一个数组的平均值*/
function arr_ave($array){
	$n=count($array);
	$sum=0;
	for($i=0;$i<$n;$i++){
	$sum+=$array[$i];
	}
	return $sum/$n;
}

/* 逐项生成一个数组的绝对值*/
function arr_abs($array){
	$out=array();
	for($i=0;$i<count($array);$i++){
	$out[$i]=abs($array[$i]);
	}
	return $out;
}

/* 向量内积*/
function arr_mul($array,$b){
	$out=array();
	for($i=0;$i<count($array);$i++){
	$out[$i]=$array[$i]*$b[$i];
	}
	return $out;
}
?>