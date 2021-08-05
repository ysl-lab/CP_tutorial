#!/bin/env bash
source ../../rerun/include.sh


dataLabels=(
" FE " 
" PE " 
" -TdS " 
)

figLabels=(
"$\Delta G\$" 
"$\Delta H\$" 
"$ -T\Delta S\$" 
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
  filter $ic > TABLE1_cluster$ic.dat
done
