#!/bin/env python

import numpy as np
import cPickle as cpk

def CalcInterEnergy (clusterID, ComplexEnergies, SoluteEnergies, SolutionEnergies):
  complex_cluster_energy = ComplexEnergies[clusterID]
  solute_cluster_energy = SoluteEnergies[clusterID]
  solution_cluster_energy = SolutionEnergies[clusterID]
#  tol = 0.01
  inter_energy = dict()
  inter_energy_std = dict()
  for key in complex_cluster_energy.keys(): 
     ene = complex_cluster_energy[key]-solute_cluster_energy[key]-solution_cluster_energy[key]
     inter_energy[key] = np.mean(ene,axis=0)             # average energy
     inter_energy_std[key] = np.std(ene,axis=0,ddof=1) 
  ene = complex_cluster_energy["PE"]-solute_cluster_energy["PE"] + solute_cluster_energy["EElr"]
  inter_energy["PE_rest"] = np.mean(ene,axis=0)
  inter_energy_std["PE_rest"] = np.std(ene,axis=0,ddof=1)
  ene = complex_cluster_energy["LJ"]-solute_cluster_energy["LJ"]
  inter_energy["LJ_rest"] = np.mean(ene,axis=0)
  inter_energy_std["LJ_rest"] = np.std(ene,axis=0,ddof=1)
  ene = complex_cluster_energy["EEsr"]-solute_cluster_energy["EEsr"]
  inter_energy["EEsr_rest"] = np.mean(ene,axis=0)
  inter_energy_std["EEsr_rest"] = np.std(ene,axis=0,ddof=1)
  ene = complex_cluster_energy["EElr"]
  inter_energy["EElr_rest"] = np.mean(ene,axis=0)
  inter_energy_std["EElr_rest"] = np.std(ene,axis=0,ddof=1)

  ene = solute_cluster_energy["EE"]-solute_cluster_energy["EElr"]
  inter_energy["EE_Pvac"] = np.mean(ene,axis=0)
  inter_energy_std["EE_Pvac"] = np.std(ene,axis=0,ddof=1)

  ene = solute_cluster_energy["PE"]-solute_cluster_energy["EElr"]
  inter_energy["PE_Pvac"] = np.mean(ene,axis=0)
  inter_energy_std["PE_Pvac"] = np.std(ene,axis=0,ddof=1)


#  assert (inter_energy['PE'] - inter_energy['PE2']) < tol  # Check if PE and PE2 are idential within some tolence
  return inter_energy, inter_energy_std
    
    
####################################################################################
import sys

complexPickleFile = sys.argv[1]
solutePickleFile = sys.argv[2]
solutionPickleFile = sys.argv[3]
#energyInfo, energyTerms = LoadEnergyConfigInfo(energyConfigFile)
complexEnergies = np.load(complexPickleFile)
soluteEnergies = np.load(solutePickleFile)
solutionEnergies = np.load(solutionPickleFile)


#Do not know why the order of the keys is different after CalcInterEnergy is called.
#print "#",
#for k in complexEnergies[1].keys():
#   print "%10s "%k,
#print 
print "#",
InterEnergy, InterEnergy_std =  CalcInterEnergy(1, complexEnergies, soluteEnergies, solutionEnergies)
for k in InterEnergy.keys():
   print "%10s "%k,
print


for ic in complexEnergies.keys():
   InterEnergy, InterEnergy_std =  CalcInterEnergy(ic, complexEnergies, soluteEnergies, solutionEnergies)
   print "# ClusterID = %d, average energy"%ic
   for k in InterEnergy.keys():
      print "%10.3f "%InterEnergy[k],
#      print "%10s %10.3f "%(k, InterEnergy[k]),
   print

print "#",
for k in InterEnergy.keys():
   print "%10s "%k,
print
for ic in complexEnergies.keys():
   InterEnergy, InterEnergy_std =  CalcInterEnergy(ic, complexEnergies, soluteEnergies, solutionEnergies)
   print "# ClusterID = %d, std"%ic
   for k in InterEnergy.keys():
      print "%10.3f "%InterEnergy_std[k],
   print

#print "#",
#for k in energyTerms:
#  print "%10s "%k,
#print 

#print " ",
#for k in energyTerms: 
#  print "%10.3f "%clusterEnergies[k],
#print 
