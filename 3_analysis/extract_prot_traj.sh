#!/bin/bash
#NUM_PROD=1
NUM_PROD=2
#NUM_PROD=3
NUM_REP=17
CURR_DIR=`pwd`
TRJCONV="gmx_mpi trjconv"
TRAJ_DIR="../2_production/"

# reference structure
TPR="../2_production/prod01/start0.tpr"

for (( ip=1; ip<=$NUM_PROD; ip++ )); do
  cd $CURR_DIR

  dir_out=`printf "prod%02d" $ip`
  if [ ! -e $dir_out ]; then
    mkdir $dir_out
  fi

  for (( ir=12; ir<$NUM_REP; ir++ )); do

    xtc_src=`printf "%s/prod%02d/prod%d.xtc" $TRAJ_DIR $ip $ir`
    #xtc_out=`printf "%s/prod%d_50_100ns.xtc" $dir_out $ir`
    xtc_out=`printf "%s/prod%d_150_200ns.xtc" $dir_out $ir`
    #xtc_out=`printf "%s/prod%d_250_300ns.xtc" $dir_out $ir`
    
    if [ -e $xtc_out ]; then
      echo "**Warning: $xtc_out already exists, skip it"
      continue
    else
      echo 3 1 | $TRJCONV -s $TPR -f $xtc_src -o $xtc_out -fit rot+trans -b 50000
    fi 
  done
done
