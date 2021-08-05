#!/bin/env bash
INI_REP=10  # neutral replica start
NUM_REP=14  # neutral replica end
PROT=cTEMPPROT
TC=trjconv_mpi
GA=g_angle_mpi

[[ ! -e s1$PROT ]] && mkdir s1$PROT
[[ ! -e s2$PROT ]] && mkdir s2$PROT

function calcDihed_s1 () {
  xtc=../../s1/3_analysis/0_100ns.xtc
  ndx=../../dPCA_3D/50_100ns/index.ndx
  out1=s1${PROT}/s1${PROT}_all.xvg
  out2=s1${PROT}/s1${PROT}_all.txt
  $GA -f $xtc -n $ndx -type dihedral -all -ov $out1
  rm -rf angdist.xvg
  tail -n +13 $out1 > $out2
  rm $out1
}

function calcDihed_s2 () {
  xtc=../../s2/3_analysis/0_100ns.xtc
  ndx=../../dPCA_3D/50_100ns/index.ndx
  out1=s2${PROT}/s2${PROT}_all.xvg
  out2=s2${PROT}/s2${PROT}_all.txt
  $GA -f $xtc -n $ndx -type dihedral -all -ov $out1 
  rm -rf angdist.xvg
  tail -n +13 $out1 > $out2
  rm $out1
}

calcDihed_s1 
calcDihed_s2
