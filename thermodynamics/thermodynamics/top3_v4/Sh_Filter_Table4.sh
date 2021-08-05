#!/bin/env bash

source ../../rerun/include.sh


dataLabels=(
#" Nonbonded " 
" LJ "
" EE "
#" Bonded " 
" Bond_P " 
" Angle_P " 
" Dih_P " 
" Imp_P " 
)

figLabels=(
"$\Delta H^{\mathrm{LJ}}\$" 
"$\Delta H^{\mathrm{EE}}\$" 
"$\Delta H_{\mathrm{P}}^{\mathrm{bond}}\$" 
"$\Delta H_{\mathrm{P}}^{\mathrm{angle}}\$" 
"$\Delta H_{\mathrm{P}}^{\mathrm{dih.}}\$" 
"$\Delta H_{\mathrm{P}}^{\mathrm{imp.}}\$" 
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
  filter $ic > TABLE4_cluster$ic.dat
done
