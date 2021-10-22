# Generate the dihedral index file
#bash VMD_GenPhiPsiIndices.sh

# Write the dpca.ndx file for PCA use 
# 16 dihedrals, each dihedral has (sin, cos) elements, totally we
# have 16x2 = 32 coordinates, we need at least 32/3 = 11 atoms to store
# these coordinates
#CHANGE prot, NT, and TIME depending on your system

NT=17
prot=cTEMPPROT
for (( it=12; it<$NT; it++ )); do
  inp=raw_traj/s1$prot/prod${it}_TIME.xtc #reordered_all.xtc
  odname=dihed_traj/s1$prot
  [[ ! -e $odname ]] && mkdir -p $odname
  out=$odname/s1${prot}_${it}.trr
  gmx_mpi angle -f $inp  -n index.ndx -or $out -type dihedral
  rm -rf angdist.xvg

  inp=raw_traj/s2$prot/prod${it}_TIME.xtc #reordered_all.xtc
  odname=dihed_traj/s2$prot
  [[ ! -e $odname ]] && mkdir -p $odname
  out=$odname/s2${prot}_${it}.trr
  gmx_mpi angle -f $inp  -n index.ndx -or $out -type dihedral
  rm -rf angdist.xvg
done


