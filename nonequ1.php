<?php
$fp=fopen("../minimize.xyz","r");
list($natom)=fscanf($fp,"%d");
$fdebug=fopen("debug.log","w");

fgets($fp);
for($i=0;$i<$natom;$i++){
	list($type[$i],$x[$i],$y[$i])=fscanf($fp,"%s%f%f");
//	fprintf($fdebug,"$type[$i]\t$x[$i]\t$y[$i]\n");
}
//neiborlist
$nC=0;
for($i=0;$i<$natom;$i++){
	if($type[$i]!="C")continue;
	$nnei[$nC]=0;
	$nnei1[$nC]=0;
	$nnei2[$nC]=0;
		$nnei3[$nC]=0;
	for($j=0;$j<$natom;$j++){
		if($type[$j]!="N")continue;
		if(within($x[$i]-$x[$j],$y[$i]-$y[$j],1.7))$nnei[$nC]++;
		else {
			if(within($x[$i]-$x[$j],$y[$i]-$y[$j],2.7))$nnei1[$nC]++;
			else 
			{
				if(within($x[$i]-$x[$j],$y[$i]-$y[$j],3))$nnei2[$nC]++;
				else {
					if(within($x[$i]-$x[$j],$y[$i]-$y[$j],4))$nnei3[$nC]++;
				}
			}
		}
	}
	$nC++;//fprintf($fdebug,"$nC\n");
}
//不等价的C个数	
$nonequ=array();
for($i=0;$i<$nC;$i++){
	$flag=1;//fprintf($fdebug,"$nnei[$i]\t");
	for($j=0;$j<count($nonequ);$j++){
		if($nnei[$i]==$nonequ[$j][0]&&$nnei1[$i]==$nonequ[$j][1]&&$nnei2[$i]==$nonequ[$j][2]&&$nnei3[$i]==$nonequ[$j][3]){$flag=0;break;}
	}
	if($flag)array_push($nonequ,array($nnei[$i],$nnei1[$i],$nnei2[$i],$nnei3[$i]));
}
echo count($nonequ);
function within($x,$y,$r){
	return $x*$x+$y*$y<$r*$r;
}
?>
