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

TC=50_100ns    # Time segment from dPCA
#TC=150_200ns    # Time segment from dPCA
#TC=250_300ns    # Time segment from dPCA
#TC=350_400ns    # Time segment from dPCA


NX=50
NY=50
NZ=50

XMIN=-6.0
XMAX=6.0
YMIN=-6.0
YMAX=6.0
ZMIN=-6.0
ZMAX=6.0

#################################################

# gmxdump eigenvec.trr file
gmxdump_mpi -f ../../dPCA_3D/${TC}/phipsi/covar/eigenvec.trr > Out_eigenvec_${prot}.out

icount=0
for d in ${dir[@]}; do
   cd ${d}
      for (( i=0; i<${end_time} ; i=i+${sep_time} )) ; do
         echo 'Processing' ${d} ${i} $((${i}+${sep_time}))
         python $pydir/Py_xvg_ranges.py ${d}_all.txt ${i} $((${i}+${sep_time}))
         python ../dPCA_project_onto_pc123.py Out_${i}_$((${i}+${sep_time}))ns.xvg ../Out_eigenvec_${prot}.out
         input=`printf "Out_%d_%dns_pc1pc2pc3.txt" ${i} $((${i}+${sep_time}))`
         den=${input/.txt/.den}
         python $pydir/Py_CalcDensit3D.py $input $XMIN $XMAX $YMIN $YMAX $ZMIN $ZMAX $NX $NY $NZ > $den    
      done
   cd ..
   icount=$((${icount}+1))
done
