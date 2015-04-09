<?
function getboxrange($path){
$file=fopen($path,'r');
for($i=0;$i<5;$i++)fscanf($file,"");
list($xlo,$xhi)=fscanf($file,"%f %f");
list($ylo,$yhi)=fscanf($file,"%f %f");
list($zlo,$zhi)=fscanf($file,"%f %f");
    return array($xlo,$xhi,$ylo,$yhi,$zlo,$zhi);
}
function getxrange($path){
        $file=fopen($path,'r');
    list($n)=fscanf($file,"%d");
    fscanf($file,"");
    $xmin=100000;$xmax=-100000;
    $ymin=100000;$ymax=-100000;
    $zmin=100000;$zmax=-100000;
    for($i=0;$i<$n;$i++){
        list($label,$x,$y,$z)=fscanf($file,"%s %f %f %f");//pay attention to the %s
        $xmin=min($x,$xmin);$xmax=max($x,$xmax);
        $ymin=min($y,$ymin);$ymax=max($y,$ymax);
        $zmin=min($z,$zmin);$zmax=max($z,$zmax);
    }
    return array($xmin,$xmax,$ymin,$ymax,$zmin,$zmax);
}


list($xlo,$xhi,$ylo,$yhi,$zlo,$zhi)=getboxrange("$fileBoxRange");
list($xlo0,$xhi0,$ylo0,$yhi0,$zlo0,$zhi0)=getxrange("$fileXRange");
if(!$xp){
	$xlo=$xlo0;$xhi=$xhi0;
}
if(!$yp){
	$ylo=$ylo0;$yhi=$yhi0;
}
if(!$zp){
	$zlo=$zlo0;$zhi=$zhi0;
}
$lx=$xhi-$xlo;$ly=$yhi-$ylo;$lz=$zhi-$zlo;
if($method=="nvt"){
	if($lx/$deta/2<$upP)exit("upP is too large!\n");
}
if($method=="muller"){
	if($lx/$deta/4<$upP)exit("upP is too large!\n");
}
$fixl1=$xlo-$deta;$fixl2=$fixl1+$deta*$wfix;
    $fixr2=$xhi+$deta;$fixr1=$fixr2-$deta*$wfix;
        if(!$hdeta)$hdeta=$deta;
    $hotl1[0]=$fixl2;$hotl2[0]=$hotl1[0]+$hdeta;
    $hotr2[0]=$fixr1;$hotr1[0]=$hotr2[0]-$hdeta;
    $cold[0]="cold0";
    $hot[0]="hot0";

    for($i=1;$i<$nstat;$i++){
        $hotl1[$i]=$hotl2[$i-1];$hotl2[$i]=$hotl1[$i]+$hdeta;
        $hotr2[$i]=$hotr1[$i-1];$hotr1[$i]=$hotr2[$i]-$hdeta;
        $cold[$i]="cold".$i;
        $hot[$i]="hot".$i;
    }
$downP=$upP;
$down11=$xlo+$deta*$downP;
$down12=$down11+$deta;
$down22=$xhi-$deta*$downP;
$down21=$down22-$deta;
$up12=$xlo+$lx/2-$deta*$upP;
$up11=$up12-$deta;
$up21=$xlo+$lx/2+$deta*$upP;
$up22=$up21+$deta;
$lp=$up11-$down11;

    	if($enforceThick)$zfactor=$lz/$thick;
    	else $zfactor=1;
    	$S=$ly*$lz;

$icoldl=$xlo;
$icoldr=$icoldl+$nswap*$deta;
$ihotl=$xlo+$lx/2;
$ihotr=$ihotl+$nswap*$deta;
?>