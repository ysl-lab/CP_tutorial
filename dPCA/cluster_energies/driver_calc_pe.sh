#!/bin/bash
CLUSTERS= #7

pop=populations.txt
assign=assignments.txt

ene1=rerun/s1/complex/energy.xvg
ene2=rerun/s1/protein/energy.xvg
ene3=rerun/s1/solvent/energy.xvg

phipsi_notreordered= #../raw_traj/s1caaAAAA/s1caaAAAA_all.xvg

#phipsi=s1cA6_reordered.xvg
#group_dir=../reordered/s1cA6
#trjcat_order=order.txt

complex=rerun/s1/complex/Out_deltaG_deltaH_deltaS.txt
protein=rerun/s1/protein/Out_deltaG_deltaH_deltaS.txt
solvent=rerun/s1/solvent/Out_deltaG_deltaH_deltaS.txt

log1=rerun/s1/complex/rerun.log
log2=rerun/s1/protein/rerun.log
log3=rerun/s1/solvent/rerun.log

pydir=/cluster/tufts/ylin12_5/Diana/Cyclic_peptides/seq_struct/python


# USE THESE FOR NON-REORDERED 

# Calculate the free energy of each cluster based on population
function calc_FE () {
   python $pydir/Py_FE_from_pop.py $pop # output is Out_FE_diff.txt
}

# Calculate phipsi angles for each cluster and make ramachandran plots
function calc_phipsi () {
   if [ ! -e phipsi ]; then
      mkdir phipsi
   fi
   python $pydir/Py_notreordered_phipsi_for_clusters.py $assign $phipsi_notreordered
   mv Out_cluster*_phipsi.txt phipsi
   sed -e "s/NS=CLUSTERS/NS=$CLUSTERS/" calc_rama.sh > tmp.sh
   bash tmp.sh
}

# Calculate average PE and the PE for all clusters
function calc_PE () {
   python $pydir/Py_notreordered_assign_energy_for_clusters.py $assign $ene1 $pop 
   mv Out_aver* Out_clust* rerun/s1/complex
   python $pydir/Py_notreordered_assign_energy_for_clusters.py $assign $ene2 $pop 
   mv Out_aver* Out_clust* rerun/s1/protein
   python $pydir/Py_notreordered_assign_energy_for_clusters.py $assign $ene3 $pop 
   mv Out_aver* Out_clust* rerun/s1/solvent
}

# Calculate deltaG deltaH deltaS for complex protein and solvent
function energies () {
   python $pydir/Py_deltaG_deltaH_clusters.py rerun/s1/complex/Out_average_PE.txt Out_FE_diff.txt 
   mv Out_delta* rerun/s1/complex
   python $pydir/Py_deltaG_deltaH_clusters.py rerun/s1/protein/Out_average_PE.txt Out_FE_diff.txt 
   mv Out_delta* rerun/s1/protein
   python $pydir/Py_deltaG_deltaH_clusters.py rerun/s1/solvent/Out_average_PE.txt Out_FE_diff.txt 
   mv Out_delta* rerun/s1/solvent
}

# Combine all energy files
function combine_energies () {
   python $pydir/Py_combine_energies_dH_dS.py $complex $protein $solvent
}

# Decomposition of deltaH
function decompose_deltaH () {
   python $pydir/Py_analyze_energies_from_log.py $log1
   mv Out_deco* rerun/s1/complex
   python $pydir/Py_analyze_energies_from_log.py $log2
   mv Out_deco* rerun/s1/protein
   python $pydir/Py_analyze_energies_from_log_sol.py $log3
   mv Out_deco* rerun/s1/solvent
}

# Decomposing energies for each cluster
function decompose_for_clusters () {
   python $pydir/Py_notreordered_decomposition_for_clusters.py $assign rerun/s1/complex/Out_decomposed_energies_from_log.txt $pop
   mv Out_average_decomposed.txt Out_cluster* rerun/s1/complex/
   python $pydir/Py_notreordered_decomposition_for_clusters.py $assign rerun/s1/protein/Out_decomposed_energies_from_log.txt $pop
   mv Out_average_decomposed.txt Out_cluster* rerun/s1/protein/
   python $pydir/Py_notreordered_decomposition_for_clusters.py $assign rerun/s1/solvent/Out_decomposed_energies_from_log.txt $pop
   mv Out_average_decomposed.txt Out_cluster* rerun/s1/solvent/
}

function top_five_clusters () {
   python $pydir/Py_dH_dS_top_five_clusters.py Out_combined_deltas.txt
}

calc_FE
calc_phipsi
calc_PE
energies
combine_energies
decompose_deltaH
decompose_for_clusters
top_five_clusters

# USE THESE FOR REORDERED
#python $pydir/Py_FE_from_pop.py $pop # output is Out_FE_diff.txt
#python $pydir/Py_reordered_indexes.py $assign $group_dir $trjcat_order
#python $pydir/Py_reordered_phipsi.py Out_mapped_indexes.txt $phipsi_notreordered
#python $pydir/Py_reordered_energies.py Out_mapped_indexes.txt $ene $pop
#python $pydir/Py_deltaG_deltaH_clusters.py Out_average_PE.txt Out_FE_diff.txt 
