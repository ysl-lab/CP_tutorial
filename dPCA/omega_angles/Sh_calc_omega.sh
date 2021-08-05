#!/bin/bash
#
# Aidan Fike
# June 27, 2019
#
# Program to obtain the trajectory of each omega angle of the cyclic peptide. 

NR=16
prodnum=$1

cp ../../../s1/1_initial/equil_npt.gro .
vmd -e VMD_GenOmegaIndex.sh #Compute the index file to indicate the omega dihedrals

#Compute the dihedrals for both s1 and s2
for ((i=0; i<=${NR} ; i++ )); do 
   g_angle_mpi -f ../../../s1/2_production/prod0${prodnum}/prod${i}.xtc -n index_omega.ndx -ov prod${i}_omega.xvg -all -type dihedral
done

mkdir s1
mv prod*xvg s1
rm angdist* \#angdist*

for ((i=0; i<=${NR} ; i++ )); do
   g_angle_mpi -f ../../../s2/2_production/prod0${prodnum}/prod${i}.xtc -n index_omega.ndx -ov prod${i}_omega.xvg -all -type dihedral
done

mkdir s2
mv prod*xvg s2
rm angdist* \#angdist*

python Py_searchForCis.py > cisOut.txt
python Py_findFirst.py > first.txt
