#!/bin/bash

PROT=cTEMPPROT
NSs1=S1NUMCLUSTERS
NSs2=S2NUMCLUSTERS

for ((i=0; i<=$NSs1; i++)); do
   python calcFreeEnergy_6mers_v2.py s1${PROT}_phipsi/cluster$i.txt s1${PROT}_cluster$i.png
done

for ((i=0; i<=$NSs2; i++)); do
   python calcFreeEnergy_6mers_v2.py s2${PROT}_phipsi/cluster$i.txt s2${PROT}_cluster$i.png
done
