#!/bin/bash

pydir=/cluster/tufts/ylin12_5/Diana/Cyclic_peptides/seq_struct/python

nclusters=
ndx=../../index_omega.ndx

for (( i=1; i<=${nclusters}; i++ )); do
   g_angle_mpi -f cluster${i}.xtc -n $ndx -ov cluster${i}.xvg -type dihedral -all
   python $pydir/Py_omega_distributions.py cluster${i}.xvg
done

rm angdist* \#angdist*
