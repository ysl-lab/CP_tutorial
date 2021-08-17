#!/bin/env bash

NRs1=S1NUMCLUSTERS
NRs2=S2NUMCLUSTERS
PROT=cTEMPPROT
TC=trjconv_mpi
GA=g_angle_mpi

function calcDihed () {
  local ip=$1
  local ir=$2
  local rep=$3

  out_dir=${rep}${PROT}_phipsi
  xtc=../cluster_traj/${rep}${PROT}/cluster${ir}.xtc
  ndx=../index.ndx
  [[ ! -e $out_dir ]] && mkdir -p $out_dir

  out=${out_dir}/cluster${ir}.xvg
  $GA -f $xtc -n $ndx -type dihedral -all -ov $out &> ${out_dir}/log_cluster${ir}.txt
  rm -rf angdist*.xvg
}

function calcDihedOld () {
  local ip=$1
  local ir=$2

  out_dir1=s1${PROT}_phipsi
  out_dir2=s2${PROT}_phipsi
  #xtc=../raw_traj/s1$PROT/all.xtc
  xtc1=../cluster_traj/s1${PROT}/cluster${ir}.xtc
  xtc2=../cluster_traj/s2${PROT}/cluster${ir}.xtc
  ndx=../index.ndx
  [[ ! -e $out_dir1 ]] && mkdir -p $out_dir1
  [[ ! -e $out_dir2 ]] && mkdir -p $out_dir2

  #out=$out_dir/s1_all.xvg
  out1=${out_dir1}/cluster${ir}.xvg
  out2=${out_dir2}/cluster${ir}.xvg
  #out=$out_dir/s1_reordered.xvg
  $GA -f $xtc1 -n $ndx -type dihedral -all -ov $out1 &> ${out_dir1}/log_cluster${ir}.txt
  $GA -f $xtc2 -n $ndx -type dihedral -all -ov $out2 &> ${out_dir2}/log_cluster${ir}.txt
  rm -rf angdist*.xvg
}


for (( ir=0; ir<=$NRs1; ir++ )); do
   calcDihed 01 $ir s1
done

for (( ir=0; ir<=$NRs2; ir++ )); do
   calcDihed 01 $ir s2
done
