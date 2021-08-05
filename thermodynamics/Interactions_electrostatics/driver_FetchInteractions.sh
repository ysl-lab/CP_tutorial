
function GetIntMat () {
  local it=$1
  local ic=$2
  local comp=$3

  inp=ENERGY_rep${it}_cluster${ic}_solute.txt
  printf "#%s\n" "NH1 CO1 NH2 CO2 NH3 CO3 NH4 CO4 NH5 CO5 NH6 CO6" 

  for g1 in NH1 CO1 NH2 CO2 NH3 CO3 NH4 CO4 NH5 CO5 NH6 CO6; do
    for g2 in NH1 CO1 NH2 CO2 NH3 CO3 NH4 CO4 NH5 CO5 NH6 CO6; do
      seq1=$comp:$g1-$g2
      seq2=$comp:$g2-$g1
      ene1=`cat $inp  | grep $seq1 `
      ene2=`cat $inp  | grep $seq2 `
      #echo "$seq1=========="
      #echo "RAW: $ene1"
      ene=0
      [[ $ene1 != "" ]] &&  ene=`echo $ene1 | awk '{print $2}'`
      #echo "FOUND: $ene"
      #echo "RAW: $ene2"
      [[ $ene2 != "" ]] &&  ene=`echo $ene2 | awk '{print $2}'`
      #echo "FOUND: $ene"
      printf "%s  " $ene 
    done
    printf "\n"
  done
}


for ((it=12; it<=16; it++ )); do
  for (( ic=8; ic<=12; ic++ )); do
    for comp in LJ-14 LJ-SR Coul-14 Coul-SR; do 
      echo "Calculating it=$it, ic=$ic, comp=$comp ..."
      out=MAT_rep${it}_cluster${ic}_$comp.txt
      GetIntMat $it $ic $comp > $out
    done
  done
done
