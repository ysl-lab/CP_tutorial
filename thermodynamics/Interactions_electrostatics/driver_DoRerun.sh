#!/bin/bash

IT=12
NT=16
IC=8
NC=12
GPP=grompp_mpi
CWD=`pwd`

function writeJobScript () {
  local xtc=$1
  local it=$2
  local ic=$3

cat << EOF > submit${it}_${ic}.job
#!/bin/bash
#SBATCH -p m4
#SBATCH -n 8 
#SBATCH --account=yushan
#SBATCH --ntasks-per-core=1
#SBATCH -J RR${it}_${ic} 
#SBATCH -o runout.%j
#SBATCH -e runerr.%j
#SBATCH -t 7-00:00:00

#module add openmpi/1.8.2
#export PATH="$PATH:/cluster/tufts/ylin12_2/HONGTAO/local/gromacs467_plumed21/bin"
#export PATH="$PATH:/cluster/tufts/ylin12_2/HONGTAO/local/gromacs467_plumed21/gromacs467/bin"
#source ~/.bash_profile

mpirun mdrun_mpi -v -rerun $xtc -s rerun.tpr -deffnm rep${it}_cluster${ic}_solute
EOF
}


function CalcClusterEnergies () {
  cd $CWD
  dname_traj=traj_cluster

#  make the tpr file
  $GPP -f files/solute.mdp -c files/solute.gro -p files/solute2.top -n files/index.ndx -o rerun.tpr &> grompp.log 

# Rerun each trajectory
  for ((it=$IT; it<=$NT; it++)); do
    for ((ic=$IC; ic<=$NC; ic++)); do
      xtc=../traj_cluster/rep${it}_cluster${ic}_solute.xtc
      writeJobScript $xtc $it $ic
      sbatch submit${it}_$ic.job
    done
  done
}


CalcClusterEnergies
