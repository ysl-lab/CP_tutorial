#!/bin/env bash
source ../../rerun/include.sh


dataLabels=(
" PE_P " 
" PE_W " 
" PE_PW " 
)

figLabels=(
"$\Delta H_{\mathrm{P}}$" 
"$\Delta H_{\mathrm{W}}$" 
"$\Delta H_{\mathrm{PW}}$" 
)

function filter () {
  local ic=$1
  inp=../top3/THERMODYNAMICS_cluster$ic.dat

  for idx in ${!dataLabels[*]}; do
    lab=${figLabels[$idx]}
    dlab=${dataLabels[$idx]}
    val=`grep  "$dlab" $inp`
    echo "\"$lab\"    $val  "
  done
}

for ((ic=$IC; ic<=$NC; ic++ )); do
  filter $ic > TABLE2_cluster$ic.dat
done
