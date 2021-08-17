PROT=
TCAT=trjcat_mpi #plumed2
TRJCONV=trjconv_mpi
G_ANGLE=g_angle_mpi
FF=rsff1
SOL=tip4pew
SIM=BE-2D-Meta


#cp ../cluster_cNPF_${FF}_${SOL}_bemeta_phipsichi18_200_300ns_Boltzmann_Reweighting/prot.gro . 


# Generate the dihedral index file
#bash VMD_GenPhiPsiIndices.sh

# Write the dpca.ndx file for PCA use 
# 16 dihedrals, each dihedral has (sin, cos) elements, totally we
# have 16x2 = 32 coordinates, we need at least 32/3 = 11 atoms to store
# these coordinates
#$TRJCONV -f prot.gro -n dpca.ndx -s prot.gro -o dpca.gro -e 0

NT=17
prot=cTEMPPROT
for (( it=12; it<$NT; it++ )); do
  inp=raw_traj/s1$prot/prod${it}_TIME.xtc #reordered_all.xtc
  odname=dihed_traj/s1$prot
  [[ ! -e $odname ]] && mkdir -p $odname
  out=$odname/s1${prot}_${it}.trr
  $G_ANGLE -f $inp  -n index.ndx -or $out -type dihedral
  rm -rf angdist.xvg

  inp=raw_traj/s2$prot/prod${it}_TIME.xtc #reordered_all.xtc
  odname=dihed_traj/s2$prot
  [[ ! -e $odname ]] && mkdir -p $odname
  out=$odname/s2${prot}_${it}.trr
  $G_ANGLE -f $inp  -n index.ndx -or $out -type dihedral
  rm -rf angdist.xvg
done


