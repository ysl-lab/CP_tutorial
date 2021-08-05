#!/bin/env bash
source ../../rerun/include.sh


dataLabels=(
" -TdSconf " 
" -TdSsolv "
)

figLabels=(
"$ -T\Delta S^{\mathrm{conf}}_{\mathrm{P}}$" 
"$ -T\Delta S_{\mathrm{W}}$"
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
  filter $ic > TABLE3_cluster$ic.dat
done
