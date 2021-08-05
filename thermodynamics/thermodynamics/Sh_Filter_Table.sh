#!/bin/env bash


dataLabels=(
" FE " 
" PE " 
" PE_P " 
" PE_W " 
" PE_PW " 
" LJ " 
" LJ_P " 
" LJ_W " 
" LJ_PW "  
" EE " 
" EE_P " 
" EE_W " 
" EE_PW "  
" Bond_P "
" Angle_P "
" Dih_P "
" Imp_P "
" -TdS " 
" -TdSconf " 
" -TdSsolv ")

figLabels=(
"$\Delta G\$" 
"$\Delta H\$" 
"$\Delta H_{P}\$" 
"$\Delta H_{W}\$" 
"$\Delta H_{PW}\$" 
"$\Delta H^{LJ}\$" 
"$\Delta H_{P}^{LJ}\$" 
"$\Delta H_{W}^{LJ}\$" 
"$\Delta H_{PW}^{LJ}\$" 
"$\Delta H^{EE}\$" 
"$\Delta H_{P}^{EE}\$" 
"$\Delta H_{W}^{EE}\$" 
"$\Delta H_{PW}^{EE}\$" 
"$\Delta H_{P}^{bond}\$" 
"$\Delta H_{P}^{angle}\$" 
"$\Delta H_{P}^{diehd}\$" 
"$\Delta H_{P}^{improper}\$" 
"-T$\Delta S\$" 
"-T$\Delta S_{conf}\$" 
"-T$\Delta S_{solv}\$"
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

for ((ic=13; ic<=17; ic++ )); do
  filter $ic > FILTERED_cluster$ic.dat
done
