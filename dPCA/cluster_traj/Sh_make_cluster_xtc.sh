#!/bin/bash

xtc=../../raw_traj/TEMPPROT/TEMPPROT_all.xtc 
gro=../../prot.gro
ndx=cluster.ndx

#clusters=( 1 2 3 4 5 ) 
max_cluster=NUMCLUSTERS
firstClust=0

##############################################

for ((i=$firstClust;i<=$max_cluster;i++)); do
   echo $i 0 | trjconv_mpi -f $xtc -s $gro -fr $ndx -o cluster${i}.xtc &> log_cluster${i}.txt
done
