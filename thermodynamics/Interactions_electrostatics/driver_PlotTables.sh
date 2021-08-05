#!/bin/env bash

for (( ic=8; ic<=12; ic++ )); do
  echo "Ploting cluster $ic ..."
  python Py_GenTable.py MAT_MEAN_cluster$ic.txt MAT_SEM_cluster$ic.txt TABLE_cluster$ic.png
  python Py_GenTable_combine.py MAT_MEAN_cluster${ic}_Coul-14.txt MAT_SEM_cluster${ic}_Coul-14.txt MAT_MEAN_cluster${ic}_LJ-14.txt MAT_SEM_cluster${ic}_LJ-14.txt TABLE_cluster${ic}_14.png
  python Py_GenTable_combine.py MAT_MEAN_cluster${ic}_Coul-SR.txt MAT_SEM_cluster${ic}_Coul-SR.txt MAT_MEAN_cluster${ic}_LJ-SR.txt MAT_SEM_cluster${ic}_LJ-SR.txt TABLE_cluster${ic}_Sr.png
  python Py_GenTable_combine.py MAT_MEAN_cluster${ic}_Coul-SR.txt MAT_SEM_cluster${ic}_Coul-SR.txt MAT_MEAN_cluster${ic}_Coul-14.txt MAT_SEM_cluster${ic}_Coul-14.txt TABLE_cluster${ic}_ALL_Coul.png
  python Py_GenTable_combine.py MAT_MEAN_cluster${ic}_LJ-SR.txt MAT_SEM_cluster${ic}_LJ-SR.txt MAT_MEAN_cluster${ic}_LJ-14.txt MAT_SEM_cluster${ic}_LJ-14.txt TABLE_cluster${ic}_ALL_LJ.png
done
