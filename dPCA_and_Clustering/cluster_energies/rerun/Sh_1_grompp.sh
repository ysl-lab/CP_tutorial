#!/bin/bash

traj=s1
traj_dir=../../../../../../${traj}/3_analysis/prod02
mdp_dir=../../mdps
top_dir=../../tops/${traj}

mkdir ${traj}
cd ${traj}

   cp ../../../../../${traj}/1_initial/equil_npt.gro .

   ### complex ###
   complex=${traj}_complex.xtc
   mkdir complex
   cd complex
      echo 'Processing' ${traj} 'complex'
      cp $traj_dir/$complex .
      echo 0 | trjconv_mpi -f ../equil_npt.gro -s ../equil_npt.gro -dump 0 -o complex.gro &> trjconv.log
      grompp_mpi -f $mdp_dir/complex.mdp -c complex.gro -p $top_dir/complex.top -o rerun.tpr &> grompp.log
      sed -e "s/TEMPXTC/$complex/" ../../submit.job > submit.job
      if [ -f rerun.tpr ] && [ -f $complex ]; then
         sbatch submit.job
      fi 
   cd ..

   ### protein ###
   protein=${traj}_prot.xtc
   mkdir protein
   cd protein
      echo 'Processing' ${traj} 'protein'
      cp $traj_dir/$protein .
      echo 1 | trjconv_mpi -f ../equil_npt.gro -s ../equil_npt.gro -dump 0 -o protein.gro &> trjconv.log
      grompp_mpi -f $mdp_dir/protein.mdp -c protein.gro -p $top_dir/protein.top -o rerun.tpr &> grompp.log
      sed -e "s/TEMPXTC/$protein/" ../../submit.job > submit.job
      if [ -f rerun.tpr ] && [ -f $protein ]; then
         sbatch submit.job
      fi
   cd ..

   ### solvent ###
   solvent=${traj}_sol.xtc
   mkdir solvent
   cd solvent
      echo 'Processing' ${traj} 'solvent'
      cp $traj_dir/$solvent .
      echo 13 | trjconv_mpi -f ../equil_npt.gro -s ../equil_npt.gro -dump 0 -o solvent.gro &> trjconv.log
      grompp_mpi -f $mdp_dir/solvent.mdp -c solvent.gro -p $top_dir/solvent.top -o rerun.tpr &> grompp.log
      sed -e "s/TEMPXTC/$solvent/" ../../submit.job > submit.job
      if [ -f rerun.tpr ] && [ -f $solvent ]; then
         sbatch submit.job
      fi
   cd ..

cd ..
