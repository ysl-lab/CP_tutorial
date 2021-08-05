#!/bin/bash

FR=12
LR=16
FC=8
LC=12

function combine () {
for ((i=$FC; i<=$LC; i++)); do
   grep "Coulomb (SR)" ENERGY_*_cluster${i}_* > COULOMB_ENERGY_cluster${i}.txt
   grep "Potential" ENERGY_*_cluster${i}_* > POTENTIAL_ENERGY_cluster${i}.txt
done
}

function calc_ave () {
for ((i=$FC; i<=$LC; i++)); do
   ave1=`awk '{ total += $2 } END { print total/NR }' POTENTIAL_ENERGY_cluster"${i}".txt`
   #ave2=`awk '{ total += $3 } END { print total/NR }' COULOMB_ENERGY_cluster"${i}".txt`
   std1=`awk '{total += $2; totalsq += $2*$2} END {print sqrt(totalsq/NR - (total/NR)**2) }' POTENTIAL_ENERGY_cluster"${i}".txt`
   #std2=`awk '{total += $3; totalsq += $3*$3} END {print sqrt(totalsq/NR - (total/NR)**2) }' COULOMB_ENERGY_cluster"${i}".txt`
   echo "cluster$i = $ave1 $std1"
   #echo "cluster$i = $ave2 $std2"
done
}

#combine
calc_ave
