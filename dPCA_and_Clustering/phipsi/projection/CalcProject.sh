#!/bin/env bash
#Change PROT depending on your system
#Calculates PC projections
PROT=cTEMPPROT
DNAME=pc1_pc2_pc3

ndx=../../dpca.ndx
gro=../../dpca.gro
vec=../covar/eigenvec.trr 

[[ ! -e $DNAME ]] && mkdir $DNAME

function project_all () {
  local dname=$1

  xtc=../../dihed_traj/all.trr
  out1=$dname/all_1.xvg
  out2=$dname/all_2.xvg
  gmx_mpi anaeig -s $gro -n $ndx -v $vec -f $xtc -2d $out1 -first 1 -last 2 &> anaeig_all1.txt
  gmx_mpi anaeig -s $gro -n $ndx -v $vec -f $xtc -2d $out2 -first 1 -last 3 &> anaeig_all2.txt

  xtc=../../dihed_traj/s1${PROT}_all.trr 
  out1=$dname/s1${PROT}_1.xvg
  out2=$dname/s1${PROT}_2.xvg
  gmx_mpi anaeig -s $gro -n $ndx -v $vec -f $xtc -2d $out1 -first 1 -last 2 &> anaeig_s11.txt
  gmx_mpi anaeig -s $gro -n $ndx -v $vec -f $xtc -2d $out2 -first 1 -last 3 &> anaeig_s12.txt

  xtc=../../dihed_traj/s2${PROT}_all.trr 
  out1=$dname/s2${PROT}_1.xvg
  out2=$dname/s2${PROT}_2.xvg
  gmx_mpi anaeig -s $gro -n $ndx -v $vec -f $xtc -2d $out1 -first 1 -last 2 &> anaeig_s21.txt
  gmx_mpi anaeig -s $gro -n $ndx -v $vec -f $xtc -2d $out2 -first 1 -last 3 &> anaeig_s22.txt
}


#This function is left here in case anyone wants to project individual replicas
function project_rep () {
  local nt=$1
  local prot=$2

  outdir=$DNAME/$prot
  [[ ! -e $outdir ]] && mkdir -p $outdir
  
  for (( it=0; it<$nt; it++ )); do
    xtc=../../dihed_traj/$prot/prod$it.trr
    out=$outdir/prod$it.xvg
    gmx_mpi anaeig -s $gro -n $ndx -v $vec -f $xtc -2d $out -first 1 -last 2 
  done
}


project_all $DNAME




