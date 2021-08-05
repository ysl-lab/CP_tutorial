#!/bin/env bash

source include.sh



function GetEnergy () {
  local atoms=$1
  cd $CWD

  [[ ! -e $atoms ]] && echo "*ERROR: can not find directory $atoms " && exit 1
  cd $atoms

  for ((it=$IT; it<=$NT; it++ )); do
    inp=rep$it.edr
    out=Energy_rep$it
    printf "%d\n" $(seq 1 100) | $GE -f $inp -o $out
  done
}


GetEnergy complex
GetEnergy solute
GetEnergy solvent 
