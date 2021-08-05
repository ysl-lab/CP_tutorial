#!/bin/env bash

source ../../rerun/include.sh


dataLabels=(
" LJ_W "
" EE_W "
)

figLabels=(
"$\Delta H_{\mathrm{W}}^{\mathrm{LJ}}\$" 
"$\Delta H_{\mathrm{W}}^{\mathrm{EE}}\$" 
)

function filter () {
  local ic=$1
  inp=THERMODYNAMICS_cluster$ic.dat
  #inp=../top3/THERMODYNAMICS_cluster$ic.dat

  for idx in ${!dataLabels[*]}; do
    lab=${figLabels[$idx]}
    dlab=${dataLabels[$idx]}
    val=`grep  "$dlab" $inp`
    echo "\"$lab\"    $val  "
  done
}

for ((ic=$IC; ic<=$NC; ic++ )); do
  filter $ic > TABLE5_cluster$ic.dat
done
