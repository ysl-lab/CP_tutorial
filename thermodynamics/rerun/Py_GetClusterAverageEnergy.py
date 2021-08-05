#!/bin/env python

import numpy as np
import cPickle as cpk

def LoadEnergyConfigInfo (inp):
   """ Read the energy configuration file for calculation. """
   lines = [line for line in open(inp,'r').readlines() if not line.startswith("#")]

   einfo = dict()
   eterms = []
   for line in lines:
     if len(line.strip()):
       k, v = line.split('=')
       einfo[k.strip()] = [ i.strip() for i in v.split('+') ]
       eterms.append(k.strip())

   return einfo, eterms


def CalcClusterEnergy (clusterID, energyInfo, allClusterEnergies):
  clusterEnergyTerms = allClusterEnergies[clusterID]
  tol = 0.01
  
  clusterEnergies = dict()
  for key, vals in energyInfo.iteritems():
    ene = 0.0
    for v in vals: ene += clusterEnergyTerms[v]
    clusterEnergies[key] = ene

  assert (clusterEnergies['PE'] - clusterEnergies['PE2']) < tol  # Check if PE and PE2 are idential within some tolence
  return clusterEnergies
    
    
####################################################################################
import sys

energyConfigFile = sys.argv[1]
clusterEnergyPickleFile = sys.argv[2]
clusterID = int(sys.argv[3])
energyInfo, energyTerms = LoadEnergyConfigInfo(energyConfigFile)
allClusterEnergies = np.load(clusterEnergyPickleFile)
clusterEnergies =  CalcClusterEnergy(clusterID, energyInfo, allClusterEnergies)

print "# ClusterID = %d, EnergyFile = %s"%(clusterID, clusterEnergyPickleFile)
print "#",
for k in energyTerms:
  print "%10s "%k,
print 

print " ",
for k in energyTerms: 
  print "%10.3f "%clusterEnergies[k],
print 
