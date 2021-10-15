#!/bin/env bash

NT=1             # number of trajectories

function gen_cluster_gmx_index () {
  dname=cluster_gmx_ndx
  [[ ! -e $dname ]] && mkdir -p $dname
  for ((it=0; it<$NT; it++ )); do
    inp=assignments.txt
    out=cluster.ndx
    python GenGromacsIndex.py $inp > $out
  done
}

gen_cluster_gmx_index



