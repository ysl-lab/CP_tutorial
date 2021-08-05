#!/bin/env bash
source ../../rerun/include.sh


dataLabels=(
" PE_P " 
" PE_W " 
" PE_PW " 
)

figLabels=(
"$\Delta H_{P}\$" 
"$\Delta H_{W}\$" 
"$\Delta H_{PW}\$" 
)

function filter () {
  local ic=$1
  inp=THERMODYNAMICS_cluster$ic.dat

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
