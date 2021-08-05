#!/bin/env bash

source ../rerun/include.sh

function writeParameterFile () {
  local ir=$1
  local ic=$2

cat << EOF > parameters
TRJ="../../traj_cluster/rep${ir}_cluster${ic}_solute.xtc"
TOP="../../tops/solute.top"
OUTDIR=output

# Use the name of the working directory as a base name for the output files ( you can change this to e. g. to name="my_project" if you prefer ) 
NAME=`pwd | awk 'BEGIN{FS="/"}{print $(NF)}'`

# Specify the number of  bins for bonds, angles and dihedrals (torsions) you want PARENT.x to use for building the 1D and 2D histograms for entropy calculation.
BBINS1D=50
ABINS1D=50
DBINS1D=50
BBINS2D=50
ABINS2D=50
DBINS2D=50

#  Specify the name of the backbone atoms as in the .top file (can be left empty, but gives more accurate results)
BACKBONE_ATOMS="CA C N H O"
EOF
}

for ((it=$IT; it<=$NT; it++ )); do
  for ((ic=$IC; ic<=$NC; ic++)); do
    cd $CWD
    dname=rep${it}_cluster$ic
    [[ -e $dname ]] && echo "**ERROR: directory $dname already exit. Terminated " && exit 1
    mkdir $dname
    cd $dname
    cp ../submit.job .
    writeParameterFile $it $ic
    sbatch submit.job
  done
done
