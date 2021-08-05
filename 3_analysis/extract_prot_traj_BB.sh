#!/bin/bash
NUM_PROD=4
NUM_REP=17
CURR_DIR=`pwd`
TRJCONV="trjconv_mpi"
TRAJ_DIR="../2_production/"


# reference structure
TPR="../2_production/prod01/start0.tpr"

#for (( ip=1; ip<=$NUM_PROD; ip++ )); do
#for (( ip=2; ip<=$NUM_PROD; ip++ )); do
#for (( ip=3; ip<=$NUM_PROD; ip++ )); do
for (( ip=4; ip<=$NUM_PROD; ip++ )); do
  cd $CURR_DIR

  dir_out=`printf "prod%02d" $ip`
  if [ ! -e $dir_out ]; then
    mkdir $dir_out
  fi

  for (( ir=12; ir<$NUM_REP; ir++ )); do

    xtc_src=`printf "%s/prod%02d/prod%d.xtc" $TRAJ_DIR $ip $ir`
    #xtc_out=`printf "%s/prod%d.xtc" $dir_out $ir`
    #xtc_out=`printf "%s/prod%d_BB_50_100ns.xtc" $dir_out $ir`
    #xtc_out=`printf "%s/prod%d_BB_150_200ns.xtc" $dir_out $ir`
    #xtc_out=`printf "%s/prod%d_BB_250_300ns.xtc" $dir_out $ir`
    xtc_out=`printf "%s/prod%d_BB_350_400ns.xtc" $dir_out $ir`
    
    if [ -e $xtc_out ]; then
      echo "**Warning: $xtc_out already exists, skip it"
      continue
    else
      #echo 3 5 | $TRJCONV -s $TPR -f $xtc_src -o $xtc_out -fit rot+trans -b 50000
      #echo 3 5 | $TRJCONV -s $TPR -f $xtc_src -o $xtc_out -fit rot+trans -b 150000
      #echo 3 5 | $TRJCONV -s $TPR -f $xtc_src -o $xtc_out -fit rot+trans -b 250000
      echo 3 5 | $TRJCONV -s $TPR -f $xtc_src -o $xtc_out -fit rot+trans -b 350000
    fi 
  done
done
