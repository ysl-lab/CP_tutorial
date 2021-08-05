#!/bin/bash

PROT=cTEMPPROT
traj=( s1 s2 )

for i in ${traj[@]}; do 
   cd ${i}${PROT}
      #trjcat_mpi -f ${i}${PROT}_10.trr ${i}${PROT}_11.trr ${i}${PROT}_12.trr ${i}${PROT}_13.trr ${i}${PROT}_14.trr -cat -o ../${i}${PROT}_all.trr &> ${i}.txt
      #trjcat_mpi -f ${i}${PROT}_11.trr ${i}${PROT}_12.trr ${i}${PROT}_13.trr ${i}${PROT}_14.trr ${i}${PROT}_15.trr -cat -o ../${i}${PROT}_all.trr
      trjcat_mpi -f ${i}${PROT}_12.trr ${i}${PROT}_13.trr ${i}${PROT}_14.trr ${i}${PROT}_15.trr ${i}${PROT}_16.trr -cat -o ../${i}${PROT}_all.trr &> ${i}.txt
      #trjcat_mpi -f ${i}${PROT}_13.trr ${i}${PROT}_14.trr ${i}${PROT}_15.trr ${i}${PROT}_16.trr ${i}${PROT}_17.trr -cat -o ../${i}${PROT}_all.trr
      #trjcat_mpi -f ${i}${PROT}_14.trr ${i}${PROT}_15.trr ${i}${PROT}_16.trr ${i}${PROT}_17.trr ${i}${PROT}_18.trr -cat -o ../${i}${PROT}_all.trr
   cd ..
done

trjcat_mpi -f s1${PROT}_all.trr s2${PROT}_all.trr -cat -o all.trr &> combined.txt
