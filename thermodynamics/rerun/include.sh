
IT=12   # Index for first neutral replica
NT=16   # Index for last neutral replica
IC=1    # Index for first cluster output
NC=48    # Index for last cluster output

module load openmpi/1.8.2
export PATH="$PATH:/cluster/tufts/ylin12_2/HONGTAO/local/gromacs467_plumed21/gromacs467/bin"
TOPDIR=../tops
MDPDIR=../tops
CWD=`pwd`
GE=g_energy_mpi
GPP=grompp_mpi
TC=trjconv_mpi
