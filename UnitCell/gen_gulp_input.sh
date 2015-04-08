#!/bin/bash
about(){
  echo "gen_gulp_input.sh v1.0 8/28/2010 Jeff Doak jwd686@u.northwestern.edu
This program reads in a gulp output file from standard input and writes a gulp
or vasp formatted input file to standard output.
Command-Line Arguments:
  -h | --help : Prints this document.
  -v | --vasp : Create a vasp-formatted input.
  -g | --gulp : Create a gulp-formatted input."
}

gulp_file(){
  #This function creates a gulp formatted input file from the results
  #of gulp calculations.
  local a=`sed -n "s/^ *a = *\([0-9]*[.][0-9]*\) *alpha = *[0-9]*[.][0-9]*/\1/p" /dev/fd/0`
  local alpha=`sed -n "s/^ *a = *[0-9]*[.][0-9]* *alpha = *\([0-9]*[.][0-9]*\)/\1/p" /dev/fd/0`
  local b=`sed -n "s/^ *b = *\([0-9]*[.][0-9]*\) *beta  = *[0-9]*[.][0-9]*/\1/p" /dev/fd/0`
  local beta=`sed -n "s/^ *b = *[0-9]*[.][0-9]* *beta  = *\([0-9]*[.][0-9]*\)/\1/p" /dev/fd/0`
  local c=`sed -n "s/^ *c = *\([0-9]*[.][0-9]*\) *gamma = *[0-9]*[.][0-9]*/\1/p" /dev/fd/0`
  local gamma=`sed -n "s/^ *c = *[0-9]*[.][0-9]* *gamma = *\([0-9]*[.][0-9]*\)/\1/p" /dev/fd/0`
  local num_atoms=`sed -n "s/^ *Total number atoms\/shells = *\([0-9][0-9]*\)/\1/p" /dev/fd/0`
  echo "cell"
  echo "$a $b $c $alpha $beta $gamma"
  #grep "Formula = " /dev/fd/0 | sed "s/ *[A-Z=][a-z]*/ /g"
  echo "Fractional"
  local flag=0; local count=0
  let local temp=$num_atoms+4
  cat /dev/fd/0 |
  while read line
  do
    if echo $line | grep -q "Final fractional coordinates of atoms :"
    then
      flag=1
    elif [ $flag -eq 1 ] && [ $count -le $temp ]
    then
      echo $line | sed -n "s/^ *[0-9][0-9]* *\([A-Z][a-z]*\) *c *\([0-9]*[.][0-9]*\) *\([0-9]*[.][0-9]*\) *\([0-9]*[.][0-9]*\) *[0-9]*[.][0-9]*/\1 \2 \3 \4/p"
      (( count++ ))
    fi
  done
  echo "output"
}

vasp_file(){
  #This function creates a vasp-formatted input file from a gulp output file
  #sent to standard input.
  local num_atoms=`sed -n "s/^ *Total number atoms\/shells = *\([0-9][0-9]*\)/\1/p" /dev/fd/0`
  echo "title"
  echo "1.0"
  local flag=0; local count=0
  let temp=3
  cat /dev/fd/0 |
  while read line
  do
    if echo $line | grep -q "Final Cartesian lattice vectors (Angstroms) :"
    then
      flag=1
    elif [ $flag -eq 1 ] && [ $count -le $temp ]
    then
      echo $line | sed -n "s/\(-*[0-9]*[.][0-9]*\) *\(-*[0-9]*[.][0-9]*\) *\(-*[0-9]*[.][0-9]*\)/\1 \2 \3/p"
      (( count++ ))
    fi
  done
  grep "Formula = " /dev/fd/0 | sed "s/ *[A-Z=][a-z]*/ /g"
  echo "Direct"
  flag=0; count=0
  let temp=$num_atoms+4
  cat /dev/fd/0 |
  while read line
  do
    if echo $line | grep -q "Final fractional coordinates of atoms :"
    then
      flag=1
    elif [ $flag -eq 1 ] && [ $count -le $temp ]
    then
      echo $line | sed -n "s/^ *[0-9][0-9]* *\([A-Z][a-z]*\) *c *\([0-9]*[.][0-9]*\) *\([0-9]*[.][0-9]*\) *\([0-9]*[.][0-9]*\) *[0-9]*[.][0-9]*/\2 \3 \4 \1/p"
      (( count++ ))
    fi
  done
}

case "$1" in
  -h|--help)
    about
    ;;
  -v|--vasp)
    vasp_file
    ;;
  -g|--gulp)
    gulp_file
    ;;
  *)
    echo "Please put either -v or -g as a command-line argument."
    ;;
esac
exit
