#!/bin/bash

dir=/cluster/tufts/ylin12/aidan/cyclic/6mersCorrectgro5/template

TOP=

cp -r $dir/02* .
cp -r $dir/03* .
cp -r $dir/04* .

cd 02_enemin
   sed -i.bak "s:TEMP:$TOP:" submit.job
   echo 'Submitting 02_enemin'
   jobid0=$(sbatch submit.job | awk '{print $4}')
cd ..

cd 03_equil
   cd res_nvt
      sed -i.bak "s:TEMP:$TOP:" submit.job
      echo 'Submitting 03_equil_nvt'
      jobid1=$(sbatch --dependency=afterok:$jobid0 submit.job | awk '{print $4}')
      jobid0=$jobid1
   cd ..
   cd res_npt
      sed -i.bak "s:TEMP:$TOP:" submit.job
      echo 'Submitting 03_equil_npt'
      jobid1=$(sbatch --dependency=afterok:$jobid0 submit.job | awk '{print $4}')
      jobid0=$jobid1
   cd ..
cd ..

cd 04_prod
   cd equil_nvt
      sed -i.bak "s:TEMP:$TOP:" submit.job
      echo 'Submitting 04_prod_nvt'
      jobid1=$(sbatch --dependency=afterok:$jobid0 submit.job | awk '{print $4}')
      jobid0=$jobid1
   cd ..
   cd equil_npt
      sed -i.bak "s:TEMP:$TOP:" submit.job
      echo 'Submitting 04_prod_npt'
      jobid1=$(sbatch --dependency=afterok:$jobid0 submit.job | awk '{print $4}')
      jobid0=$jobid1
   cd ..
cd ..


