#!/bin/env bash

source ../rerun/include.sh


dataLabels=(
" EE " 
" EE14_P " 
" EEsr_P " 
" EElr_P " 
" EEsr_W " 
" EElr_W " 
" EEsr_PW "  
" EElr_PW "  
)

figLabels=(
"$\Delta H_{EE}\$" 
"$\Delta H_{P}^{EE(1-4)}\$" 
"$\Delta H_{P}^{EE(SR)}\$" 
"$\Delta H_{P}^{EE(LR)}\$" 
"$\Delta H_{W}^{EE(SR)}\$" 
"$\Delta H_{W}^{EE(LR)}\$" 
"$\Delta H_{PW}^{EE(SR)}\$" 
"$\Delta H_{PW}^{EE(LR)}\$" 

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
