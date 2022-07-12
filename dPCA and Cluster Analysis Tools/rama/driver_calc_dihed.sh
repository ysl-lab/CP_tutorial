#!/bin/env bash
#Change PROT below to match your system
#NR determines number of clusters plotted
NR=5
PROT=cGNSRV


function calcDihed () {
  local ip=$1
  local ir=$2

  out_dir1=s1${PROT}_phipsi
  out_dir2=s2${PROT}_phipsi
  
  xtc1=../cluster_traj/s1${PROT}/cluster${ir}.xtc
  xtc2=../cluster_traj/s2${PROT}/cluster${ir}.xtc
  ndx=../index.ndx
  [[ ! -e $out_dir1 ]] && mkdir -p $out_dir1
  [[ ! -e $out_dir2 ]] && mkdir -p $out_dir2


  out1=${out_dir1}/cluster${ir}.xvg
  out2=${out_dir2}/cluster${ir}.xvg

  gmx_mpi angle -f $xtc1 -n $ndx -type dihedral -all -ov $out1 &> ${out_dir1}/log_cluster${ir}.txt 
  gmx_mpi angle -f $xtc2 -n $ndx -type dihedral -all -ov $out2 &> ${out_dir2}/log_cluster${ir}.txt 
  rm -rf angdist*.xvg
} 



for (( ir=1; ir<=$NR; ir++ )); do
   calcDihed 01 $ir
done
