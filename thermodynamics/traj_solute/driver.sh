NUM_PROD=1
INI_REP=12
NUM_REP=16
BEGIN=50000      # ps 
CWD=`pwd`
TC="trjconv_mpi"

TRAJ_DIR="../../../../s1/2_production"
TPR="../../../../s1/2_production/prod01/start0.tpr"
GRO="../../../../s1/2_production/prod01/prod0.gro"


[[ ! -e "solute.gro" ]] && echo 1 | $TC -s $TPR -f $GRO -o solute.gro


for (( ip=1; ip<=$NUM_PROD; ip++ )); do
  cd $CWD

  dir_out=`printf "prod%02d" $ip`
  [[ ! -e $dir_out ]] && mkdir $dir_out

  for (( ir=$INI_REP; ir<=$NUM_REP; ir++ )); do
    xtc_src=`printf "%s/prod%02d/prod%d.xtc" $TRAJ_DIR $ip $ir`
    xtc_out=`printf "%s/rep%d.xtc" $dir_out $ir`
    
    if [ -e $xtc_out ]; then
      echo "**Warning: $xtc_out already exists, skip it"
      continue
    else
      echo 1 | $TC -s $TPR -f $xtc_src -o $xtc_out -b $BEGIN -pbc mol
    fi 
  done

done
