#!/bin/bash

NREP=

function gen_bemeta_temp () {
  output=$1

cat BE_2D_MetaD_Input.dat << EOF > $1 



EOF

}


function driver () {
  gen_bemeta_temp tmp.dat
  for (( ir=0; ir<$NREP; ir++ )); do
    out=`printf "bemeta.dat.%d" $ir`
    cvx=cv$((ir*2))
    cvy=cv$((ir*2+1))
    sed -e "s/cvX/$cvx/g; s/cvY/$cvy/g; " tmp.dat > $out
  done
  rm tmp.dat
}

driver
