#!/bin/env bash

source ../rerun/include.sh

################ Free Energy from Population ###################
function calcFE() {
  for ((ic=$IC; ic<=$NC; ic++)); do
    out=FE_cluster$ic.dat
    fes="FE      "
    TdS="-TdSconf"
    for ((it=$IT; it<=$NT; it++)); do
      echo "Process cluster $ic, replica $it ..."
      inp1=../traj_cluster/population_rep$it.txt 
      inp2=../entropy/ENTROPY_rep$it.txt
      [[ ! -e $inp1 ]] && echo "**ERROR: can not file file $inp1. Terminated. " && exit 1
      [[ ! -e $inp2 ]] && echo "**ERROR: can not file file $inp2. Terminated. " && exit 1
      fe=`python Py_CalcFEFromPop.py $inp1 $ic`
      mTdS=`python Py_CalcEntropy.py $inp2 $ic`
      fes="$fes    $fe"
      TdS="$TdS    $mTdS"
    done
    echo $fes > $out
    echo $TdS >> $out
  done
}


################ Enthalpy ###################
function calcEnthalpy() {
#  REF="../rerun/ENTHALPY_cluster_$NC.dat"
  for ((ic=$IC; ic<=$NC; ic++)); do 
     inp1="../rerun/ENTHALPY_5rep_Ave.dat"
     inp2="../rerun/ENTHALPY_5rep_Std.dat"
     out=ENTHALPY_cluster_${ic}_adjusted.dat
     python Py_GetClusterEnthalpy_adjusted.py  $inp1 $inp2 $ic > $out 
  done
}

################ All ###################
function calcThermodynamics() {
  for ((ic=$IC; ic<=$NC; ic++)); do
    echo "Calculating thermodynamics for cluster $ic"
    python Py_CalcClusterThermodynamics_adjusted.py ENTHALPY_cluster_${ic}_adjusted.dat FE_cluster$ic.dat  > THERMODYNAMICS_cluster${ic}.dat
  done
}

#calcFE
calcEnthalpy  # The enthalpy terms for output is defined in Enthalpy.info
calcThermodynamics


#bash Sh_Filter_Table1.sh
#bash Sh_Filter_Table2.sh
#bash Sh_Filter_Table3.sh
#bash Sh_Filter_Table4.sh
#python Py_GenTable.py TABLE1
#python Py_GenTable.py TABLE2
#python Py_GenTable.py TABLE3
#python Py_GenTable.py TABLE4
