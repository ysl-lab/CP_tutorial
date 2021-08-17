#!/bin/bash

NS=CLUSTERS

for ((i=1; i<=$NS; i++)); do
   python calcFreeEnergy.py phipsi/Out_cluster${i}_phipsi.txt phipsi/Out_cluster${i}_phipsi.png
done
