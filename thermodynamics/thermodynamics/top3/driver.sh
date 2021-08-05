#!/bin/env bash

source ../../rerun/include.sh

################ Free Energy from Population ###################
function calcFE() {
  for ((ic=$IC; ic<=$NC; ic++)); do
    out=FE_cluster$ic.dat
    fes="FE      "
    TdS="-TdSconf"
    for ((it=$IT; it<=$NT; it++)); do
      echo "Process cluster $ic, replica $it ..."
      inp1=../../traj_cluster/population_rep$it.txt 
      inp2=../../entropy/ENTROPY_rep$it.txt
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
  REF="../../rerun/ENTHALPY_cluster_$NC.dat"
  for ((ic=$IC; ic<=$NC; ic++)); do
     inp="../../rerun/ENTHALPY_cluster_$ic.dat"
     out=ENTHALPY_cluster$ic.dat
     python Py_GetClusterEnthalpy.py  $inp $REF > $out 
  done
}

################ All ###################
function calcThermodynamics() {
  for ((ic=$IC; ic<=$NC; ic++)); do
    echo "Calculating thermodynamics for cluster $ic"
    python Py_CalcClusterThermodynamics.py ENTHALPY_cluster$ic.dat FE_cluster$ic.dat  > THERMODYNAMICS_cluster$ic.dat
  done
}

#calcFE
#calcEnthalpy  # The enthalpy terms for output is defined in Enthalpy.info
#calcThermodynamics


bash Sh_Filter_Table1.sh # Table1 = G H S
bash Sh_Filter_Table2.sh # Table2 = H breakdown
bash Sh_Filter_Table3.sh # Table3 = S breakdown
bash Sh_Filter_Table4.sh # Table4 = bonded-nonbonded
bash Sh_Filter_Table5.sh # Table5 = coulomb
bash Sh_Filter_Table6.sh # Table6 = LJ

python Py_GenTable.py TABLE1
python Py_GenTable.py TABLE2
python Py_GenTable2.py TABLE3
python Py_GenTable3.py TABLE4
python Py_GenTable3.py TABLE5
python Py_GenTable3.py TABLE6
