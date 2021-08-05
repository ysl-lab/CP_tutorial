#!/bin/bash 

SIM=1

NUM_REP=17                # number of replicas
NUM_CPU_PER_REP=8          # number cpus for each replica
REP_EX_RATE=2500           # replica exchange every 2500 steps
TEMPERATURE=300.0          # Simulation temperature
PRESSURE=1.0               # Simulation pressure
EXTEND_TIME=100000         # Extended simulation time 100 ns

NODE="m4c[]"
JOB_NAME="cTEMPPROT"
MDRUN="gmx_mpi mdrun"
GROMPP="gmx_mpi grompp"
TPBCONV="gmx_mpi tpbconv" 
TOP_DIR="../../1_initial"
TOP="cTEMPPROT_rsff2_tip3p.top"
NUM_CPU=$(($NUM_REP * $NUM_CPU_PER_REP))


#******************************************************************************#
#                                                                              #
# Generate MDP from a template file                                            #
#                                                                              #
#******************************************************************************#
function write_mdp ()
{
  local mdp_template=$1
  local output=$2

  if [ ! -e $output ]; then
    sed -e "s#TEMP#$TEMPERATURE#g; s#PRES#$PRESSURE#g" $mdp_template > $output
  else
    echo ">> File $output exists, skip the generation"
  fi
}


#******************************************************************************#
#                                                                              #
# Generate the LSF job script                                                  #
#                                                                              #
#******************************************************************************#
function write_job_script () {
  local mdrun_command=$1

  echo "#!/bin/bash"
  echo "#SBATCH -p m4"
  echo "#SBATCH -n $NUM_REP"
  echo "#SBATCH --cpus-per-task=$NUM_CPU_PER_REP"
  echo "#SBATCH --exclude=pcomp[01-44]"
  echo "#SBATCH --qos=yushan"
  echo "#SBATCH --account=yushan"
  echo "#SBATCH --ntasks-per-core=1"
  echo "#SBATCH -J $JOB_NAME"
  echo "#SBATCH -o runout.%j"
  echo "#SBATCH -e runerr.%j"
  echo "#SBATCH -t 7-00:00:00"
  echo "#SBATCH --export=ALL"
  echo " "
  echo "$mdrun_command"
}


#******************************************************************************#
#                                                                              #
# Generate TPR from GRO/MPD/TOP                                                #
#                                                                              #
#******************************************************************************#
function generate_tpr ()
{
  local gro=$1
  local mdp=$2
  local top=$3
  local tpr=$4
  local log=${5:-"grompp_log.txt"}

  if [ ! -e $tpr ]; then
    $GROMPP -v -f $mdp -p $top -c $gro -o $tpr >& $log
  else
    echo "**File $tpr exists, skip the generation"
  fi
}


#******************************************************************************#
#                                                                              #
# Check a directory, if it exists, terminate current process; otherwise        #
# create the directory                                                         #
#                                                                              #
#******************************************************************************#
function check_dir ()
{
  dir=$1
  if [ ! -e $dir ]; then
    mkdir -p $dir
  else
    echo "Directory $dir exists, skip it"
    exit 1
  fi
}



################################################################################
#                                                                              #
#                   Here starts the main process                               #
#                                                                              #
################################################################################




# NVT equilibration

if [ "$SIM" == "1" ]; then   

  mdp_template="../prod_npt_bemeta.mdp"
  prod_dir="prod01"
  gro="../../1_initial/equil_npt.gro"

  check_dir "$prod_dir"
  cd $prod_dir

  cp ../bemeta.dat* .

  cp $TOP_DIR/$TOP . 
  write_mdp $mdp_template "prod_npt.mdp" 
  write_job_script "mpiexec -np $NUM_REP $MDRUN -v -ntomp $NUM_CPU_PER_REP -plumed bemeta -multi $NUM_REP -replex $REP_EX_RATE -s start -deffnm prod" > submit.job 

  for (( irep=0; irep<$NUM_REP; irep++ )); do
    output_tpr=`printf "start%d.tpr" $irep`
    generate_tpr $gro "prod_npt.mdp" $TOP $output_tpr
  done  
else
  bash restart.sh

  prev_dir=`printf "../prod%02d" $(($SIM-1)) `
  prod_dir=`printf "prod%02d" $SIM `

  check_dir $prod_dir
  if [ $SIM -gt 9 ]; then
	echo "Problem with script, HILLS not copied"
  fi
  cd $prod_dir
  cp ../prod0$(($SIM - 1))/HILLS.* .
  cp ../prod0$(($SIM - 1))/COLVAR.* .

  cp ../bemeta.dat* .
  write_job_script "mpiexec -np $NUM_REP $MDRUN -v -ntomp $NUM_CPU_PER_REP -plumed bemeta -multi $NUM_REP -replex $REP_EX_RATE -cpi $prev_dir/prod.cpt -s start.tpr -deffnm prod" > submit.job

  for (( irep=0; irep<$NUM_REP; irep++ )); do
    output_tpr=`printf "start%d.tpr" $irep`
    tpbconv_log=`printf "tpbconv_log%d.txt" $irep`
    $TPBCONV -s $prev_dir/$output_tpr -extend $EXTEND_TIME -o $output_tpr   >& $tpbconv_log
  done  
fi

#bsub < submit.job 


