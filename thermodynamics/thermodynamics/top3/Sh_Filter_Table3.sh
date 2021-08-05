#!/bin/env bash
source ../../rerun/include.sh


dataLabels=(
" -TdSconf " 
" -TdSsolv "
)

figLabels=(
"$ -T\Delta S^{conf}_P\$" 
"$ -T\Delta S_{W}\$"
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
  filter $ic > TABLE3_cluster$ic.dat
done
