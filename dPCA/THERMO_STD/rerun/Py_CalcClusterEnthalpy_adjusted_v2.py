#!/bin/env python

import numpy as np
import cPickle as cpk

def LoadEnergies (inp):
   """ Read the energy configuration file for calculation. """

   energies = np.loadtxt(inp, dtype=float, comments="#")


   # Find the keys for builidng up a dict
   for line in open(inp,'r').readlines():
     if line.startswith("#"):
       if "PE2" in line: 
         energyLabels = line.split()[1:]
         break

   allEnergies = dict()
   ii = 0
   n = len(energyLabels)
   for ene in energies:
     assert n == len(ene)
     repEnergies = dict()
     for i in range(n): repEnergies[energyLabels[i]] = ene[i]
     
     allEnergies[ii] = repEnergies
     ii += 1

   return energyLabels, allEnergies
      
###############################################################################
import sys

complexEnergyFile = sys.argv[1]
soluteEnergyFile = sys.argv[2]
solventEnergyFile = sys.argv[3]
interEnergyFile = sys.argv[4]
complexEnergyLabels, complexEnergies = LoadEnergies (complexEnergyFile)
soluteEnergyLabels, soluteEnergies = LoadEnergies (soluteEnergyFile)
solventEnergyLabels, solventEnergies = LoadEnergies (solventEnergyFile)
interEnergyLabels, interEnergies = LoadEnergies (interEnergyFile)

nterm = len(complexEnergyLabels)
trajKeys = complexEnergies.keys()
for ii in range(nterm):
  label = complexEnergyLabels[ii]

# Complex thermodynamics
  print "%15s"%label,
  for k in trajKeys:   # Loop through all trajectories
    print "%15.3f"%complexEnergies[k][label],
  print 

# Solute thermodynamics
  if label in soluteEnergies[0]:
    print "%13s_P"%label,
    for k in trajKeys:   # Loop through all trajectories
      print "%15.3f"%soluteEnergies[k][label],
    print 

# Solvent thermodynamics
  if label in solventEnergies[0]:
    print "%13s_W"%label,
    for k in trajKeys:   # Loop through all trajectories
      print "%15.3f"%solventEnergies[k][label],
    print 
  
# Solute/Solvent thermodynamics
  if label in interEnergies[0]:
    print "%12s_PW"%label,
    for k in trajKeys:   # Loop through all trajectories
      print "%15.3f"%interEnergies[k][label],
    print

for label in interEnergies[0]:
  if "rest" in label or "Pvac" in label:
    print "%15s"%label,
    for k in trajKeys:   # Loop through all trajectories
      print "%15.3f"%interEnergies[k][label],
    print


#  if label in soluteEnergies[0] and label in solventEnergies[0]:
#    print "%12s_PW"%label,
#    for k in trajKeys:   # Loop through all trajectories
#      print "%15.3f"%(complexEnergies[k][label] - soluteEnergies[k][label] - solventEnergies[k][label]),
#    print 
 
