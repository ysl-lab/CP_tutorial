
IT=12
NT=16
IC=8
NC=12

for ((it=$IT; it<=$NT; it++)); do
  for ((ic=$IC; ic<=$NC; ic++)); do
    fname=rep${it}_cluster${ic}_solute
    echo "Processing $fname ..."
    printf "%d\n" $(seq 1 800) | g_energy_mpi -f $fname.edr -o $fname.xvg > ENERGY_$fname.txt
  done
done
