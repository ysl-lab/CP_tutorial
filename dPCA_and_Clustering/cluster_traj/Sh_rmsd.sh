#!/bin/bash

pydir=/cluster/tufts/ylin12_5/Diana/Cyclic_peptides/seq_struct/python
NMR=/cluster/tufts/ylin12_5/Diana/Cyclic_peptides/seq_struct/6mers/16_NMe15_aAAAAA/NMR_1.5_CHA.pdb
#NMR=/cluster/tufts/ylin12_5/Diana/Cyclic_peptides/seq_struct/6mers/15_NMe16_aAAAAA/NMR_1.6_CHA.pdb

nclusters=

make_ndx_mpi -f $NMR -o index.ndx << EOF
a C | a CA | a N | a O 
keep 10
q
EOF

for (( i=1; i<=${nclusters}; i++ )); do
   g_rms_mpi -s $NMR -f cluster${i}.xtc -o cluster${i}_rmsd.xvg -fit rot+trans -pbc -n index.ndx
   python $pydir/Py_average_rmsd.py cluster${i}_rmsd.xvg
done

