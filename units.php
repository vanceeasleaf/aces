<?php
$u="real";
    $boltz[$u]= 0.0019872067;
    $hplanck[$u]= 95.306976368;
    $mvv2e[$u]= 48.88821291 * 48.88821291;
    $ftm2v[$u]= 1.0 / 48.88821291 / 48.88821291;
    $mv2d[$u]= 1.0 / 0.602214179;
    $nktv2p[$u]= 68568.415;
    $qqr2e[$u]= 332.06371;
    $qe2f[$u]= 23.060549;
    $vxmu2f[$u]= 1.4393264316e4;
    $xxt2kmu[$u]= 0.1;
    $e_mass[$u]= 1.0/1836.1527556560675;
    $hhmrr2e[$u]= 0.0957018663603261;
    $mvh2r[$u]= 1.5339009481951;
    $angstrom[$u]= 1.0;
    $femtosecond[$u]= 1.0;
    $qelectron[$u]= 1.0;
    $kelvin[$u]=1;
    $timestep= 1.0;
    $ar_mass[$u]=39.95;
$u="metal";
    $boltz[$u]= 8.617343e-5;
    $hplanck[$u]= 4.135667403e-3;
    $mvv2e[$u]= 1.0364269e-4;
    $ftm2v[$u]= 1.0 / 1.0364269e-4;
    $mv2d[$u]= 1.0 / 0.602214179;
    $nktv2p[$u]= 1.6021765e6;
    $qqr2e[$u]= 14.399645;
    $qe2f[$u]= 1.0;
    $vxmu2f[$u]= 0.6241509647;
    $xxt2kmu[$u]= 1.0e-4;
    $e_mass[$u]= 0.0;    // not yet set
    $hhmrr2e[$u]= 0.0;
    $mvh2r[$u]= 0.0;
    $angstrom[$u]= 1.0;
    $femtosecond[$u]= 1.0e-3;
    $qelectron[$u]= 1.0;
    $kelvin[$u]=1;
    $timestep= 0.001;
    $ar_mass[$u]=39.95;
$u="si";
    $boltz[$u]= 1.3806504e-23;
    $hplanck[$u]= 6.62606896e-34;
    $mvv2e[$u]= 1.0;
    $ftm2v[$u]= 1.0;
    $mv2d[$u]= 1.0;
    $nktv2p[$u]= 1.0;
    $qqr2e[$u]= 8.9876e9;
    $qe2f[$u]= 1.0;
    $vxmu2f[$u]= 1.0;
    $xxt2kmu[$u]= 1.0;
    $e_mass[$u]= 0.0;    // not yet set
    $hhmrr2e[$u]= 0.0;
    $mvh2r[$u]= 0.0;
    $angstrom[$u]= 1.0e-10;
    $femtosecond[$u]= 1.0e-15;
    $qelectron[$u]= 1.6021765e-19;
    $kelvin[$u]=1;
    $timestep= 1.0e-8;
      $ar_sigma[$u]=3.405e-10;//in si
  $ar_epsilon[$u]=1.67e-21;
  $ar_mass[$u]=6.633e-26;
$u="cgs";
    $boltz[$u]= 1.3806504e-16;
    $hplanck[$u]= 6.62606896e-27;
    $mvv2e[$u]= 1.0;
    $ftm2v[$u]= 1.0;
    $mv2d[$u]= 1.0;
    $nktv2p[$u]= 1.0;
    $qqr2e[$u]= 1.0;
    $qe2f[$u]= 1.0;
    $vxmu2f[$u]= 1.0;
    $xxt2kmu[$u]= 1.0;
    $e_mass[$u]= 0.0;    // not yet set
    $hhmrr2e[$u]= 0.0;
    $mvh2r[$u]= 0.0;
    $angstrom[$u]= 1.0e-8;
    $femtosecond[$u]= 1.0e-15;
    $qelectron[$u]= 4.8032044e-10;

    $kelvin[$u]=1;
    $timestep= 1.0e-8;
  $ar_mass[$u]=6.633e-23;
$u="electron";
    $boltz[$u]= 3.16681534e-6;
    $hplanck[$u]= 0.1519829846;
    $mvv2e[$u]= 1.06657236;
    $ftm2v[$u]= 0.937582899;
    $mv2d[$u]= 1.0;
    $nktv2p[$u]= 2.94210108e13;
    $qqr2e[$u]= 1.0;
    $qe2f[$u]= 1.94469051e-10;
    $vxmu2f[$u]= 3.39893149e1;
    $xxt2kmu[$u]= 3.13796367e-2;
    $e_mass[$u]= 0.0;    // not yet set
    $hhmrr2e[$u]= 0.0;
    $mvh2r[$u]= 0.0;
    $angstrom[$u]= 1.88972612;
    $femtosecond[$u]= 0.0241888428;
    $qelectron[$u]= 1.0;
    $kelvin[$u]=1;
    $timestep= 0.001;
    $ar_mass[$u]=39.95e-3;
$u="lj";
    $boltz[$u]= 1.0;// J/K
    $hplanck[$u]= 0.18292026;  // using LJ parameters for argon
    $mvv2e[$u]= 1.0;
    $ftm2v[$u]= 1.0;
    $mv2d[$u]= 1.0;
    $nktv2p[$u]= 1.0;
    $qqr2e[$u]= 1.0;
    $qe2f[$u]= 1.0;
    $vxmu2f[$u]= 1.0;
    $xxt2kmu[$u]= 1.0;
    $e_mass[$u]= 0.0;    // not yet set
    $hhmrr2e[$u]= 0.0;
    $mvh2r[$u]= 0.0;
    $ar_mass[$u]=1.0;
    $angstrom[$u]= $angstrom["si"]/$ar_sigma["si"];
    $femtosecond[$u]= $femtosecond["si"]*sqrt($ar_epsilon["si"]/$ar_mass["si"]/$ar_sigma["si"]/$ar_sigma["si"]);
    $qelectron[$u]= 1.0;
    $kelvin[$u]=$kelvin["si"]*($boltz["si"]/$ar_epsilon["si"]);
    $timestep= 0.005;
function Unit($srcUnits,$type){
global $units;
global $angstrom;
global     $femtosecond;
global    $boltz;
global $ar_mass;
global     $kelvin;
switch($type){
case "l":return $angstrom[$units]/$angstrom[$srcUnits];break;
case "t":return $femtosecond[$units]/$femtosecond[$srcUnits];break;
case "T":return $kelvin[$units]/$kelvin[$srcUnits];break;
case "M":return $ar_mass[$units]/$ar_mass[$srcUnits];break;
case "E":return ($boltz[$units]*$kelvin[$units])/($boltz[$srcUnits]*$kelvin[$srcUnits]);break;
default: exit("unknown units type");
}
}

?>
