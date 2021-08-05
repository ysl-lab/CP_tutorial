#!/bin/bash

source include.sh
CWD=`pwd`

function writeJobScript () {
  local xtc=$1
  local it=$2

cat << EOF > submit$it.job
#!/bin/bash
#SBATCH -p m4
#SBATCH -n 2
#SBATCH --account=yushan
#SBATCH --qos=yushan
#SBATCH --ntasks-per-core=1
#SBATCH -J RR$it 
#SBATCH -o runout.%j
#SBATCH -e runerr.%j
#SBATCH -t 7-00:00:00

module add openmpi/1.8.2
export PATH="$PATH:/cluster/tufts/ylin12_2/HONGTAO/local/gromacs467_plumed21/bin"
export PATH="$PATH:/cluster/tufts/ylin12_2/HONGTAO/local/gromacs467_plumed21/gromacs467/bin"
source ~/.bash_profile
mpirun mdrun_mpi -v -rerun $xtc -s rerun.tpr -deffnm rep$it
EOF
}


function CalcEnergyRerun () {
  local atoms=$1
  cd $CWD

  dname_traj=../traj_$atoms
  dname_out=$atoms
  [[ ! -e $dname_out ]] && mkdir -p $dname_out

  cd $CWD  
  cd $dname_out

#  make the tpr file
  cp ../$dname_traj/$atoms.gro .
  cp ../$TOPDIR/$atoms.top .
  cp ../$TOPDIR/$atoms.mdp .
  $GPP -f $atoms -c $atoms -p $atoms -o rerun.tpr &> grompp.log 

 
# Rerun each trajectory
  for ((it=$IT; it<=$NT; it++)); do
    xtc=../$dname_traj/prod01/rep$it.xtc
    writeJobScript $xtc $it
    sbatch submit$it.job
  done
}

CalcEnergyRerun solute
CalcEnergyRerun complex
CalcEnergyRerun solvent
