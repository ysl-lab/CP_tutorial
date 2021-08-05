num_cv=17
for (( i=0; i<$num_cv; i++ )); do
  inp="bemeta.dat.$i"
  out="bemeta.dat.$i"
  sed -i -e '1iRESTART\'  $inp
done
