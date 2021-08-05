#!/bin/env bash

#GNUPLOT=$HOME/hongtao/local/bin/gnuplot
GNUPLOT=gnuplot
TC=trjconv
NT=1             # number of trajectories
NS=8             # number of states
FF=amber03
SOL=tip3p
STRUCT=S1


function assign () {
  dname=traj_pca
  for ((it=0; it<$NT; it++ )); do
    inp=$dname/prod$it.txt
    out=$dname/asn_prod$it.txt
    ./assign $inp $out
  done
}


function gen_cluster_gmx_index () {
  dname=cluster_gmx_ndx
  [[ ! -e $dname ]] && mkdir -p $dname
  for ((it=0; it<$NT; it++ )); do
    inp=assignments.txt
    out=cluster.ndx
    python GenGromacsIndex.py $inp > $out
  done
}

function find_cluster_phipsi () {
  for ((it=0; it<$NT; it++ )); do
    asn=traj_pca/asn_prod$it.txt
    out=traj/dihed_prod$it
    echo " >> $asn ..."
    dihed=$out.xvg
    python FindClusterPhiPsi.py $asn $dihed $out
  done
}

function concat () {
  for ((is=1; is<$NS; is++ )); do
    out=phipsi_cluster$is.txt
    rm -rf $out
    for ((it=0; it<$NT; it++ )); do
      inp=traj/dihed_prod${it}_cluster_${is}.txt
      cat $inp >> $out
    done
  done
}

function make_plot () {
  for (( is=1; is<=$NS; is++ )); do
    inp=phipsi_cluster$is.txt
    $GNUPLOT -e "INPUT='$inp'; NS='$NS'; IS='$is';" make_plots_phipsi.gplt
    convert -density 300 phipsi.eps state_$is.png
  done
}


# Clean the comments from the trajectory PCA files
#assign
#plot_traj_cluster
gen_cluster_gmx_index
#find_cluster_phipsi
#concat
#make_plot 


