#!/bin/env bash

source ../../rerun/include.sh


dataLabels=(
#" LJ " 
" LJ14_P " 
" LJsr_P " 
" LJlr_P " 
" LJsr_W " 
" LJlr_W " 
" LJsr_PW "  
" LJlr_PW "  
)

figLabels=(
#"$\Delta H_{LJ}\$" 
"$\Delta H_{P}^{LJ(1-4)}\$" 
"$\Delta H_{P}^{LJ(SR)}\$" 
"$\Delta H_{P}^{LJ(LR)}\$" 
"$\Delta H_{W}^{LJ(SR)}\$" 
"$\Delta H_{W}^{LJ(LR)}\$" 
"$\Delta H_{PW}^{LJ(SR)}\$" 
"$\Delta H_{PW}^{LJ(LR)}\$" 

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
  filter $ic > TABLE6_cluster$ic.dat
done
