#!/bin/bash
#Change PROT below to match your system
#NS decides number of states to be plotted
PROT=cTEMPPROT
NS=5

for ((i=1; i<=$NS; i++)); do
   echo 'Calculating cluster' ${i}
   python calcFreeEnergy_v2_5mers.py s1${PROT}_phipsi/cluster$i.txt s1${PROT}_cluster$i.png
   python calcFreeEnergy_v2_5mers.py s2${PROT}_phipsi/cluster$i.txt s2${PROT}_cluster$i.png
done
