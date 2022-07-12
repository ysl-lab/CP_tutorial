#!/bin/bash
#Change PROT below to match your system
#Change clusters and max_cluster below to match your number of clusters
xtc=../../raw_traj/s1cPROT/s1cPROT_all.xtc 
gro=../../prot.gro
ndx=cluster.ndx

clusters=( 1 2 3 4 5 6 7 ) 
max_cluster=7

##############################################

for c in ${clusters[@]}; do
   #echo $c,$max_cluster
   echo $max_cluster 0 | gmx_mpi trjconv -f $xtc -s $gro -fr $ndx -o cluster${c}.xtc &> log_cluster${c}.txt
   max_cluster=$((${max_cluster}-1))
done
