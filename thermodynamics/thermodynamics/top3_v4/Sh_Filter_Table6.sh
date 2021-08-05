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
#"$\Delta H_{\mathrm{LJ}}\$" 
"$\Delta H_{\mathrm{P}}^{\mathrm{LJ(1-4)}}\$" 
"$\Delta H_{\mathrm{P}}^{\mathrm{LJ(SR)}}\$" 
"$\Delta H_{\mathrm{P}}^{\mathrm{LJ(LR)}}\$" 
"$\Delta H_{\mathrm{W}}^{\mathrm{LJ(SR)}}\$" 
"$\Delta H_{\mathrm{W}}^{\mathrm{LJ(LR)}}\$" 
"$\Delta H_{\mathrm{PW}}^{\mathrm{LJ(SR)}}\$" 
"$\Delta H_{\mathrm{PW}}^{\mathrm{LJ(LR)}}\$" 

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
  filter $ic > TABLE6_cluster$ic.dat
done
