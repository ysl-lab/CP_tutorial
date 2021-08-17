#!/bin/env bash
PROT=cTEMPPROT
G_ANAEIG=g_anaeig_mpi
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
  $G_ANAEIG -s $gro -n $ndx -v $vec -f $xtc -2d $out1 -first 1 -last 2 &> anaeig_all1.txt
  $G_ANAEIG -s $gro -n $ndx -v $vec -f $xtc -2d $out2 -first 1 -last 3 &> anaeig_all2.txt

  xtc=../../dihed_traj/s1${PROT}_all.trr 
  out1=$dname/s1${PROT}_1.xvg
  out2=$dname/s1${PROT}_2.xvg
  $G_ANAEIG -s $gro -n $ndx -v $vec -f $xtc -2d $out1 -first 1 -last 2 &> anaeig_s11.txt
  $G_ANAEIG -s $gro -n $ndx -v $vec -f $xtc -2d $out2 -first 1 -last 3 &> anaeig_s12.txt

  xtc=../../dihed_traj/s2${PROT}_all.trr 
  out1=$dname/s2${PROT}_1.xvg
  out2=$dname/s2${PROT}_2.xvg
  $G_ANAEIG -s $gro -n $ndx -v $vec -f $xtc -2d $out1 -first 1 -last 2 &> anaeig_s21.txt
  $G_ANAEIG -s $gro -n $ndx -v $vec -f $xtc -2d $out2 -first 1 -last 3 &> anaeig_s22.txt

#  xtc=../../dihed_traj/s3$PROT/s3${PROT}_0.trr
#  out1=$dname/s3${PROT}_1.xvg
#  out2=$dname/s3${PROT}_2.xvg
#  $G_ANAEIG -s $gro -n $ndx -v $vec -f $xtc -2d $out1 -first 1 -last 2
#  $G_ANAEIG -s $gro -n $ndx -v $vec -f $xtc -2d $out2 -first 1 -last 3

#  xtc=../../dihed_traj/s4$PROT/s4${PROT}_0.trr
#  out1=$dname/s4${PROT}_1.xvg
#  out2=$dname/s4${PROT}_2.xvg
#  $G_ANAEIG -s $gro -n $ndx -v $vec -f $xtc -2d $out1 -first 1 -last 2
#  $G_ANAEIG -s $gro -n $ndx -v $vec -f $xtc -2d $out2 -first 1 -last 3

}



function project_rep () {
  local nt=$1
  local prot=$2

  outdir=$DNAME/$prot
  [[ ! -e $outdir ]] && mkdir -p $outdir
  
  for (( it=0; it<$nt; it++ )); do
    xtc=../../dihed_traj/$prot/prod$it.trr
    out=$outdir/prod$it.xvg
    $G_ANAEIG -s $gro -n $ndx -v $vec -f $xtc -2d $out -first 1 -last 2 
  done
}

# Projection of NMR structures
function project_nmr () {
  xtc=../../../dPCA_NMR/cNPF_NMR_dihed.trr 
  out=$DNAME/NMR.xvg
  $G_ANAEIG -s $gro -n $ndx -v $vec -f $xtc -2d $out -first 1 -last 2 
}


project_all $DNAME
#project_rep 18 s1cNPF
#project_rep 18 s2cNPF
#project_nmr







