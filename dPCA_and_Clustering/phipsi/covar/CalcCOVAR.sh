#!/bin/env bash

PROT=AA                         # protein name
NCV=1                           # number of collective variables

function calc_pca_covar () {
  xtc=../../dihed_traj/all.trr 
  gro=../../dpca.gro
  #[[ ! -e $outdir ]] && mkdir -p $outdir

  gmx_mpi covar -n ../../dpca.ndx -f $xtc -s $gro -ascii -xpm -xpma -nofit -nomwa -noref -nopbc &> log_covar.txt
}

calc_pca_covar
