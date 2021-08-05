#!/bin/env bash

source include.sh
nl=$((NC*2+3))

function GetClusterEnergy () {
  local atoms=$1
  cd $CWD
  [[ ! -e $atoms ]] && echo "*ERROR: can not find directory $atoms " && exit 1
  energyOutput="ENERGY_5rep_${atoms}.dat"
  output=$atoms/Pickle_Save_Energy.dat
 
  energyFile1=$atoms/Energy_rep12.xvg
  assignmentFile1=../assignments/assignment_rep12.txt
  energyFile2=$atoms/Energy_rep13.xvg
  assignmentFile2=../assignments/assignment_rep13.txt
  energyFile3=$atoms/Energy_rep14.xvg
  assignmentFile3=../assignments/assignment_rep14.txt
  energyFile4=$atoms/Energy_rep15.xvg
  assignmentFile4=../assignments/assignment_rep15.txt
  energyFile5=$atoms/Energy_rep16.xvg
  assignmentFile5=../assignments/assignment_rep16.txt
  energyConfigFile=EnergyInfo.$atoms
#  python Py_GetClusterEnergy_with_error_5rep_v2.py $energyFile1 $assignmentFile1 $energyFile2 $assignmentFile2 $energyFile3 $assignmentFile3 $energyFile4 $assignmentFile4 $energyFile5 $assignmentFile5  $energyConfigFile > $energyOutput
  python Py_GetClusterEnergy_with_error_5rep_v3.py $energyFile1 $assignmentFile1 $energyFile2 $assignmentFile2 $energyFile3 $assignmentFile3 $energyFile4 $assignmentFile4 $energyFile5 $assignmentFile5  $energyConfigFile $output > $energyOutput

head -n $nl ENERGY_5rep_${atoms}.dat > ENERGY_5rep_${atoms}_Ave.dat
tail -n $nl ENERGY_5rep_${atoms}.dat > ENERGY_5rep_${atoms}_Std.dat


}

GetClusterEnergy complex
GetClusterEnergy solute
GetClusterEnergy solvent

python Py_GetCluster_PW_int_Energy_v2.py complex/Pickle_Save_Energy.dat solute/Pickle_Save_Energy.dat solvent/Pickle_Save_Energy.dat > ENERGY_5rep_PW.dat
head -n $nl ENERGY_5rep_PW.dat > ENERGY_5rep_PW_Ave.dat
tail -n $nl ENERGY_5rep_PW.dat > ENERGY_5rep_PW_Std.dat
python Py_CalcClusterEnthalpy_adjusted_v2.py ENERGY_5rep_complex_Ave.dat ENERGY_5rep_solute_Ave.dat ENERGY_5rep_solvent_Ave.dat ENERGY_5rep_PW_Ave.dat > ENTHALPY_5rep_Ave.dat
python Py_CalcClusterEnthalpy_adjusted_v2.py ENERGY_5rep_complex_Std.dat ENERGY_5rep_solute_Std.dat ENERGY_5rep_solvent_Std.dat ENERGY_5rep_PW_Std.dat > ENTHALPY_5rep_Std.dat
