#!/bin/env bash

source ../../rerun/include.sh


################ All ###################
function calcThermodynamics() {
  for ((ic=$IC; ic<=$NC; ic++)); do
    echo "Calculating thermodynamics for cluster $ic"
    python Py_CalcClusterThermodynamics.py ../top3/ENTHALPY_cluster$ic.dat ../top3/FE_cluster$ic.dat  > THERMODYNAMICS_cluster$ic.dat
  done
}

calcThermodynamics


bash Sh_Filter_Table1.sh # Table1 = G H S
bash Sh_Filter_Table2.sh # Table2 = H breakdown
bash Sh_Filter_Table3.sh # Table3 = S breakdown
bash Sh_Filter_Table4.sh # Table4 = Protein breakdown (LJ, EE, bonded)
bash Sh_Filter_Table5.sh # Table5 = Water breakdown (LJ, EE)
bash Sh_Filter_Table6.sh # Table5 = Protein-Water breakdown (LJ, EE)

python Py_GenTable.py TABLE1
python Py_GenTable.py TABLE2
python Py_GenTable2.py TABLE3
python Py_GenTable3.py TABLE4
python Py_GenTable2.py TABLE5
python Py_GenTable2.py TABLE6

rm -rf Figs_v3
mkdir -p Figs_v3
mv *.png Figs_v3
