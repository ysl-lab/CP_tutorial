kb=0.008314 #kJ/mol/K
T=300.0     # Temperature
source ../rerun/include.sh


function calcEntropy () {
  local ir=$1
  local ic=$2

  NAME=rep${ir}_cluster${ic}
  H=`tail -n 1 $NAME/output/${NAME}_MIST.txt | awk '{print $5}'`
  S=$(echo "$kb*$H" | bc -l)
  TS=$(echo "$T*$S" | bc -l)

  printf "%8d %8.3f %8.3f\n" $ic $S $TS
}

for ((it=$IT; it<=$NT; it++ )); do
  out=ENTROPY_rep$it.txt
  for ((ic=$IC; ic<=$NC; ic++)); do
    echo "Collecting entropy for trajectory $it, cluster $ic ..."
    [[ $ic == $IC ]] && rm -r $out
    calcEntropy $it $ic >> $out
  done
done
