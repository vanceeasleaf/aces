<?php
require_once("discript.php");
require_once("postMini.php");
if(!$langevin)$langevin=0;
if(!$nvt)$nvt=0;
@passthru("python $home/input.py $units $xp $yp $zp $dumpRate $timestep $method $kb $nktv '$masses' '$potential' $T $seed $dtime $equTime $langevin $nvt $aveRate $deta $jprofile $dumpRate $corRate $computeTc $lx $ly $lz $fourierTc $tcfactor $zfactor $gstart $jcf  $nswap $excRate $lp $S $excNum $swapEnergyRate $dumpxyz $dumpv $runTime>input1");
/* 空格*/
$s='  ';
?>
#settings
units            <?echo $units?> 
dimension        3
#newton            on

boundary       <?echo $xp?"p":"s" ?>  <?echo $yp?"p":"s" ?>   <?echo $zp?"p":"s" ?> 
atom_style        atomic
read_restart   <?echo $fileRestartMini?> 
change_box all	boundary       <?echo $xp?"p":"s" ?>  <?echo $yp?"p":"s" ?>   <?echo $zp?"p":"s" ?> 
lattice            fcc    5 #needed to define the regions
thermo <?echo $dumpRate?> 
thermo_modify     lost warn
timestep          <?echo $timestep?>


#regions and groups
<?	if($method=="nvt"){?>
region            stayl    block  <?echo $fixl1,$s,$fixl2?> INF INF INF  INF units box
<?
    for($i=0;$i<$nstat;$i++){
    echo
        "region            $cold[$i]    block   $hotl1[$i]  $hotl2[$i] INF INF INF  INF units box
group            $cold[$i]    region $cold[$i]
region            $hot[$i]    block  $hotr1[$i]   $hotr2[$i] INF INF INF  INF units box
group            $hot[$i]    region $hot[$i]
";}?>
region            re_nve    block  <?echo $hotl2[$nstat-1],$s,$hotr1[$nstat-1]?> INF INF INF  INF units box
region            stayr    block  <?echo $fixr1,$s,$fixr2?> INF INF INF  INF units box
region            stay    union  2 stayl stayr
group            stay    region stay
group            g_nve    region re_nve
region          main    union <?echo $nstat*2+1?>  re_nve <?for($i=0;$i<$nstat;$i++){echo "$cold[$i] $hot[$i] ";}?> 
group   main region main
	<?}?>
	<?if($method=="muller"||$method=="inject"){?>
region            up1    block   <? echo $up11,$s,$up12?>  INF INF INF INF  units box
region            up2    block   <? echo $up21,$s,$up22?>  INF INF INF INF  units box
region            up  union 2 up1 up2
region            down1  block  <? echo $down11,$s, $down12?>   INF INF INF INF units box
region            down2  block  <? echo $down21,$s, $down22?>   INF INF INF INF units box
region            down union 2 down1 down2
	<?}?>
region hot  block   <? echo $ihotl,$s,$ihotr?>  INF INF INF INF  units box
region cold block   <? echo $icoldl,$s,$icoldr?>  INF INF INF INF  units box
#computes
compute           ke  all  ke/atom
compute           pe  all  pe/atom
compute         stress all stress/atom virial
compute jflux all heat/flux ke pe stress
	<?if($method=="muller"||$method=="inject"){?>
compute           up_temp    all  temp/region up
compute           down_temp  all  temp/region down
		<?}?>
	<? if($method!="greenkubo"){?>
#variables
variable   jx atom (c_pe+c_ke)*vx-(c_stress[1]*vx+c_stress[4]*vy+c_stress[5]*vz)/<?echo $nktv2p[$units]?> 
variable   jy atom (c_pe+c_ke)*vy-(c_stress[4]*vx+c_stress[2]*vy+c_stress[6]*vz)/<?echo $nktv2p[$units]?> 
variable   jz atom (c_pe+c_ke)*vz-(c_stress[5]*vx+c_stress[6]*vy+c_stress[3]*vz)/<?echo $nktv2p[$units]?> 
variable          temp atom  c_ke/(1.5*<?echo $boltz[$units] ?>)
	variable te atom c_ke+c_pe
	variable jcx atom v_te*vx
	variable jcy atom v_te*vy
	variable jcz atom v_te*vz
variable          delta_temp   equal  c_up_temp-c_down_temp
	<?}?>
#init atoms to T
<?	echo $masses;?> 
<? echo $potential;?> 
reset_timestep 0
velocity  all  create  <?echo $T?> <?echo $seed?> mom yes  rot yes dist gaussian
<?
if($method=="nvt"){
    echo"
velocity stay set 0 0 0
fix getEqu  main  nvt temp $T $T $dtime
";
    }else{
    echo"
fix getEqu  all  nvt temp $T $T $dtime
";
}
?>
run <?echo $equTime?> 
unfix getEqu
reset_timestep 0
<?
if($method=="nvt"){
	if(!$langevin){
	    echo"
	fix     nve  g_nve  nve
	";
	}else{
	echo"	fix     nve  main  nve
			";
	}
    }else{
    	if($nvt){//a key to elimilate the heat of numeric
    	    echo"
fix getEqu  all  nvt temp $T $T $dtime
";
    	}else{
     echo"
fix   nve  all  nve
";
   	 }
}
?>
fix    flux_out  all  ave/time  1  <?echo $aveRate?>  <?echo $aveRate?>  c_jflux[1]  c_jflux[2] c_jflux[3] file  <?echo $fileFlux?> 
<?
if($method=="nvt"){
	if($langevin){
			$tmstat="langevin";
		$r="$seed tally yes";
	}else{
		$tmstat="nvt temp";
					$r="";
	}
for($i=0;$i<$nstat;$i++){
    echo"fix   $hot[$i]  $hot[$i] $tmstat  $Thi $Thi  $dtime $r
fix   $cold[$i] $cold[$i] $tmstat  $Tlo $Tlo  $dtime $r
";
if(!$langevin){
echo"
compute my$hot[$i] $hot[$i] temp/com
fix_modify $hot[$i] temp my$hot[$i]
compute my$cold[$i] $cold[$i] temp
compute_modify my$cold[$i] extra 0
fix_modify $cold[$i] temp my$cold[$i]
";
}
else{
echo"
compute my$hot[$i] $hot[$i] temp
fix_modify $hot[$i] temp my$hot[$i]
compute my$cold[$i] $cold[$i] temp
fix_modify $cold[$i] temp my$cold[$i]
";
}
}?>

variable hot equal 0<?for($i=0;$i<$nstat;$i++){echo "+f_$hot[$i]";}?> 
variable cold equal 0<?for($i=0;$i<$nstat;$i++){echo "+f_$cold[$i]";}?> 
fix     j_hot all ave/time 1 <?echo $aveRate?> <?echo $aveRate?> v_hot  v_cold file <?echo $fileNvtWork?> 
<?}?>
<? if($method!="greenkubo"){
    echo "
	fix	temp_profile    all    ave/spatial  1  $aveRate  $aveRate  x  lower  $deta      v_temp  v_jx file  $fileTempProfile  norm sample units box
	";
	/* 输出热流空间分布，不计算热导率*/
	if($jprofile){
	echo "
	dump jprofile all custom $dumpRate $fileJProfile id v_jx v_jy v_jz v_temp v_jcx v_jcy v_jcz vx vy vz x y z
		dump_modify  jprofile sort id
	";
	}
}else {

$v=$lx*$ly*$lz;
$kb=$boltz[$units];
$factor=$corRate*$timestep/($v*$kb*$T*$T)*$zfactor*$tcfactor;

/* 用傅里叶变换热流来计算热流关联函数*/
if($fourierTc){
echo "
	fix               j_out  all  ave/time  1  1  1  c_jflux[1] c_jflux[2] c_jflux[3]  file  jin.txt 
";
}

/* 用分子模拟论坛的compute扩展来计算热流关联函数，比lammps自带的更精确*/
if($computeTc){
	$rfactor=$tcfactor*$zfactor;
	if(!$gstart)$gstart=20000;
	echo "
variable          factor_ac equal 1.0
variable          factor_tc equal $rfactor
compute           tc all tc c_thermo_temp c_jflux v_factor_ac v_factor_tc x first $gstart 0 500000
fix               tc_out  all  ave/time  1  1  1  c_tc   file  $fileKappa
		";
}
if(!$fourierTc&&!$computeTc)
{
    	echo "
    		fix ss all ave/correlate $corRate $corNum $aveRate c_jflux[1] c_jflux[2] c_jflux[3] type auto ave running
variable k11 equal trap(f_ss[3])*$factor
variable k22 equal trap(f_ss[4])*$factor
variable k33 equal trap(f_ss[5])*$factor
fix output all ave/time 1  1 $aveRate v_k11  v_k22  v_k33 file $fileKappa
#fix output1 all ave/time $aveRate 1 $aveRate v_k11 v_k22 v_k33 file kappa1.txt ave running
";

/* 定时输出热流自关联函数*/
if($jcf){
	echo "
		fix out all ave/time $aveRate 1 $aveRate f_ss[3] f_ss[4] f_ss[5] mode vector file jcf*.txt
		";
		}
		}
}
if($method=="muller"){
	$Nswapbins=2*floor($lx/(2*$nswap*$deta));
	$factor=$lz/($excRate*$timestep)/(2*$S)/(1/$lp)*$zfactor*$tcfactor;
    echo "
fix delta_out  all  ave/time  1  $aveRate  $aveRate  v_delta_temp   file  $fileTempDiff
fix heat_swap   all  thermal/conductivity  $excRate  x   $Nswapbins
fix e_exchange  all  ave/time  $excRate  $excNum  $aveRate  f_heat_swap  file  $fileSwap
variable thermal_conductivity equal f_e_exchange/(1e-10+f_delta_out)*$factor
fix thermal_conductivity_out  all  ave/time  $aveRate  1   $aveRate  v_thermal_conductivity   file  $fileTherCon
";
}
if($method=="inject"){
		$factor=1/(2*$S)/(1/$lp)*$zfactor*$tcfactor;
	echo"
		fix delta_out  all  ave/time  1  $aveRate  $aveRate  v_delta_temp   file  $fileTempDiff
fix               hot   all  heat  $excRate   $swapEnergyRate   region hot
fix               cold  all  heat  $excRate  -$swapEnergyRate   region cold
variable          thermal_conductivity equal $swapEnergyRate/(1e-10+f_delta_out)*$factor
fix               thermal_conductivity_out  all  ave/time  $aveRate  1   $aveRate  v_thermal_conductivity   file  $fileTherCon
";
}

/* 定时输出dump文件并按id排序*/
if($dumpxyz){
    echo "
dump dump1 all atom $dumpRate $fileDump
dump_modify  dump1 sort id
";}

/* 定时输出速度文件用于计算速度关联函数*/
if($dumpv){
		    echo "
dump dump2 all custom $dumpRate $fileDumpVelocity type vx vy vz
				dump_modify  dump2 sort id
";
}
?>
run           <? echo $runTime?>



