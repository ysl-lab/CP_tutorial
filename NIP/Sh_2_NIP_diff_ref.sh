#!/bin/bash

prot=cTEMPPROT

runs=( s1 s2 )
dir=( s1${prot} s2${prot} )
sep_time=25    # ns
end_time=100   # ns
#sep_time=50    # ns
#end_time=200   # ns
#sep_time=50    # ns
#end_time=300   # ns


pydir=/cluster/tufts/ylin12_5/Diana/Cyclic_peptides/seq_struct/python/NIP

#################################################

## first make reference density as the average of the 4 methods (not MD) ##
st=$((${end_time}-${sep_time}))
#ref_den=Out_ave_${st}_${end_time}ns.den

#python $pydir/Py_ave_density.py s1${prot}/Out_${st}_${end_time}ns_pc1pc2pc3.den s2${prot}/Out_${st}_${end_time}ns_pc1pc2pc3.den $ref_den
#exit

## then calculate the NIP to the reference density
for d in ${dir[@]}; do
   if [ $d == 's1'${prot} ]; then 
      ref_den='s2'${prot}/Out_${st}_${end_time}ns_pc1pc2pc3.den
   else ref_den='s1'${prot}/Out_${st}_${end_time}ns_pc1pc2pc3.den
   fi
   echo $d,$ref_den
   cd ${d}
      for (( i=0; i<${end_time} ; i=i+${sep_time} )) ; do
         echo 'Processing' ${d} ${i} $((${i}+${sep_time}))
         input=`printf "Out_%d_%dns_pc1pc2pc3.den" ${i} $((${i}+${sep_time}))`
         python $pydir/Py_calcNIP_3D.py ../$ref_den $input > NIP_${i}_$((${i}+${sep_time}))ns.txt
      done
      icount=1
      for (( i=0; i<${end_time} ; i=i+${sep_time} )) ; do
         file=NIP_${i}_$((${i}+${sep_time}))ns.txt
         cont=`cat $file`
         echo '"'${i}'-'$((${i}+${sep_time}))'ns''"' $icount $cont >> NIP_${d}.txt
         icount=$((${icount}+1))
      done
   cd ..
done
