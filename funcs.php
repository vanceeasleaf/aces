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


?>