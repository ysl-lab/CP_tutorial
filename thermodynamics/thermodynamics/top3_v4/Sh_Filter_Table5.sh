#!/bin/env bash

source ../../rerun/include.sh


dataLabels=(
#" EE " 
" EE14_P " 
" EEsr_P " 
" EElr_P " 
" EEsr_W " 
" EElr_W " 
" EEsr_PW "  
" EElr_PW "  
)

figLabels=(
#"$\Delta H_{\mathrm{EE}}\$" 
"$\Delta H_{\mathrm{P}}^{\mathrm{EE(1-4)}}\$" 
"$\Delta H_{\mathrm{P}}^{\mathrm{EE(SR)}}\$" 
"$\Delta H_{\mathrm{P}}^{\mathrm{EE(LR)}}\$" 
"$\Delta H_{\mathrm{W}}^{\mathrm{EE(SR)}}\$" 
"$\Delta H_{\mathrm{W}}^{\mathrm{EE(LR)}}\$" 
"$\Delta H_{\mathrm{PW}}^{\mathrm{EE(SR)}}\$" 
"$\Delta H_{\mathrm{PW}}^{\mathrm{EE(LR)}}\$" 
)

function filter () {
  local ic=$1
  inp=THERMODYNAMICS_cluster$ic.dat
  inp=../top3/THERMODYNAMICS_cluster$ic.dat

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
