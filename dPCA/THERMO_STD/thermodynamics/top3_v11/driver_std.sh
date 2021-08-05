#!/bin/env bash

source ../../rerun/include.sh


################ All ###################
function calcThermodynamics() {
  for ((ic=$IC; ic<=$NC; ic++)); do
    echo "Calculating thermodynamics for cluster $ic"
    python Py_CalcClusterThermodynamics_adjusted_v2.py ../ENTHALPY_cluster_${ic}_adjusted.dat ../FE_cluster$ic.dat  > THERMODYNAMICS_cluster$ic.dat
  done
}

calcThermodynamics


bash Sh_Filter_Table1.sh # Table1 = G H S
bash Sh_Filter_Table2.sh # Table2 = H breakdown
bash Sh_Filter_Table3.sh # Table3 = S breakdown
bash Sh_Filter_Table4_adjusted.sh # Table4 = Protein breakdown (LJ, EE, bonded)
bash Sh_Filter_Table5.sh # Table5 = Water breakdown (LJ, EE)

echo "$NC"
python Py_GenTable_adjusted.py TABLE1
python Py_GenTable2_adjusted.py TABLE2
python Py_GenTable2_adjusted.py TABLE3
python Py_GenTable3_adjusted.py TABLE4
python Py_GenTable_adjusted.py TABLE5

rm -rf Figs_v10
mkdir -p Figs_v10
mv *.png Figs_v10

