#!/bin/env bash

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

rm -rf Figs
mkdir -p Figs
mv *.png Figs
