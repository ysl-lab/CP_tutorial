#!/bin/env bash

source ../../rerun/include.sh


dataLabels=(
" LJ_PW "
" EE_PW "
)

figLabels=(
"$\Delta H_{\mathrm{PW}}^{\mathrm{LJ}}\$" 
"$\Delta H_{\mathrm{PW}}^{\mathrm{EE}}\$" 
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
  filter $ic > TABLE6_cluster$ic.dat
done
