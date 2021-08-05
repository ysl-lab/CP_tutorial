#!/bin/env bash

PROT=AA                         # protein name
NCV=1                           # number of collective variables
TRJCAT=trjcat_mpi #plumed2           # trjcat binary program
TRJCONV=trjconv_mpi #plumed2         # trjconv binary program
G_RMSF=g_rmsf_mpi #plumed2           # g_rmsf binary program
G_COVAR=g_covar_mpi         
G_ANAEIG=g_anaeig_mpi #plumed
XPM2PS=xpm2ps_mpi #plumed2
CWD=`pwd`


function calc_pca_covar () {
  xtc=../../dihed_traj/all.trr 
  gro=../../dpca.gro
  #[[ ! -e $outdir ]] && mkdir -p $outdir

  $G_COVAR -n ../../dpca.ndx -f $xtc -s $gro -ascii -xpm -xpma -nofit -nomwa -noref -nopbc &> log_covar.txt
}

calc_pca_covar
