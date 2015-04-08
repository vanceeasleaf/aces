<?php
include("structure.php");
$a=array();
$b=array(1,2);
$a=array_merge($a,$b);
print_r($a);
$unitCell=getUnitCell();
print count($unitCell);