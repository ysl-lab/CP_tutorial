#!/bin/env bash

source ../rerun/include.sh

function GenClusterGMXFrameIndex () {

  for ((it=$IT; it<=$NT; it++ )); do
    inp=../assignments/assignment_rep$it.txt
    pop=population_rep$it.txt
    rm -rf $pop
    for ((ic=$IC; ic<=$NC; ic++)); do
      out=index_rep${it}_cluster$ic.ndx
      echo "Processing replica $it, cluster $ic ..."
      python Py_GenGromacsIndex.py $inp $out $ic >> $pop
    done
  done
}


function WriteClusterXTC () {
  local atoms=$1

  for ((it=$IT; it<=$NT; it++ )); do
    for ((ic=$IC; ic<=$NC; ic++)); do
      xtc=../traj_$atoms/prod02/rep$it.xtc
      out=rep${it}_cluster${ic}_$atoms.xtc
      ndx=index_rep${it}_cluster${ic}.ndx
      echo "Processing trajectory $atoms, replica $it, cluster $ic ..."
      $TC -f $xtc -o $out -fr $ndx


      #echo 1 0 |$TC -f $xtc -o temp.xtc -fr $ndx -center -pbc mol -s ../rerun/$atoms/rerun.tpr
      #echo 5 0 | $TC -f temp.xtc -s ../rerun/$atoms/rerun.tpr -o $out -fit rot+trans
    done
  done

}

GenClusterGMXFrameIndex
WriteClusterXTC solute
WriteClusterXTC complex
