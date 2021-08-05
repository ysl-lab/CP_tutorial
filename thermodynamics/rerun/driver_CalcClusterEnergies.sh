#!/bin/env bash

source include.sh



function GetClusterEnergy () {
  local atoms=$1
  cd $CWD

  [[ ! -e $atoms ]] && echo "*ERROR: can not find directory $atoms " && exit 1

  for ((it=$IT; it<=$NT; it++ )); do
    energyFile=$atoms/Energy_rep$it.xvg
    assignmentFile=../assignments/assignment_rep$it.txt
    output=$atoms/Pickle_Cluster_Energy_rep$it.dat
    python Py_GetClusterEnergy.py $energyFile $assignmentFile $output

    
    for ((ic=$IC; ic<=$NC; ic++ )); do
      energyOutput="ENERGY_${atoms}_cluster_$ic.dat"
      [[ $it == $IT ]] && rm -r $energyOutput 
      python Py_GetClusterAverageEnergy.py EnergyInfo.$atoms $atoms/Pickle_Cluster_Energy_rep$it.dat $ic >> $energyOutput
    done
  done
}


GetClusterEnergy complex
GetClusterEnergy solute
GetClusterEnergy solvent 
