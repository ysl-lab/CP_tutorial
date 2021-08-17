#!/bin/bash

PROT=cTEMPPROT

function combine(){
   local fname=$1
   inp1=pc1_pc2_pc3/${fname}_1.xvg
   inp2=pc1_pc2_pc3/${fname}_2.xvg

   tail -n +13 $inp1 > tmp1 
   tail -n +13 $inp2 > tmp2

   echo "   Processing $inp1 $inp2 ..."
   python Py_combine.py tmp1 tmp2 > pc1_pc2_pc3/$fname.txt
}

combine all
combine s1${PROT}
combine s2${PROT}

rm tmp*
