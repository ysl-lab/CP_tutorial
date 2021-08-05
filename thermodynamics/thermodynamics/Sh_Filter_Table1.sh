#!/bin/env bash
source ../rerun/include.sh


dataLabels=(
" FE " 
" PE " 
" -TdS " 
" PE_P " 
" PE_W " 
" PE_PW " 
" -TdSconf " 
" -TdSsolv "
)

figLabels=(
"$\Delta G\$" 
"$\Delta H\$" 
"$ -T\Delta S\$" 
"$\Delta H_{P}\$" 
"$\Delta H_{W}\$" 
"$\Delta H_{PW}\$" 
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
  filter $ic > TABLE1_cluster$ic.dat
done
