#!/bin/env bash

source ../rerun/include.sh


dataLabels=(
" Nonbonded " 
" LJ "
" EE "
" Bonded " 
" Bond_P " 
" Angle_P " 
" Dih_P " 
" Imp_P " 
)

figLabels=(
"$\Delta H^{nonbonded}\$" 
"$\Delta H^{LJ}\$" 
"$\Delta H^{EE}\$" 
"$\Delta H_{P}^{bonded}\$" 
"$\Delta H_{P}^{bond}\$" 
"$\Delta H_{P}^{angle}\$" 
"$\Delta H_{P}^{dih}\$" 
"$\Delta H_{P}^{imp}\$" 
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
  filter $ic > TABLE4_cluster$ic.dat
done
