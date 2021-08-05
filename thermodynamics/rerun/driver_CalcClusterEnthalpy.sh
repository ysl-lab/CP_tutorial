#!/bin/env bash

source include.sh


for ((ic=$IC; ic<=$NC; ic++ )); do
  echo "Calculating thermodynamics for cluster $ic ..."
  python Py_CalcClusterEnthalpy.py ENERGY_complex_cluster_$ic.dat ENERGY_solute_cluster_$ic.dat ENERGY_solvent_cluster_$ic.dat > ENTHALPY_cluster_$ic.dat
done
