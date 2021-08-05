#!/bin/env python


import numpy as np
import sys
import cPickle as cpk 

def GetLegend (inp):
  lines = open(inp,'r').readlines()

  labels = []
  for line in lines:
    if line.startswith("@ s"):
      xx = line.split('"')
      labels.append(xx[-2]) 
  return labels


def GetEnergies (inp):
  lines = [line.split() for line in open(inp,'r').readlines() if not (line.startswith("#") or line.startswith("@"))]
  energies = np.array(lines,dtype=float)
  return energies 


#################################################################################################################
energyFile = sys.argv[1]
assignmentFile = sys.argv[2]
output = sys.argv[3]
print ">>Input Energy File: %s"%energyFile
print ">>Input Assignment File: %s"%assignmentFile
print ">>Output File: %s"%output

labels   = GetLegend(energyFile)
energies = GetEnergies(energyFile)
cids = np.loadtxt(assignmentFile, dtype=int, usecols=(4,))
assert len(cids) == len(energies)
cluster_energies = dict()

nc = np.max(cids)  # number of clusters
for ic in range(nc+1):
   u =  energies[np.where(cids==ic)]  # energy of each cluster
   aveu = np.mean(u,axis=0)           # average energy
   nf = len(u)                        # number of frames in each cluster
   print "Number of frames for cluster %d: %d"%(ic, nf)
   xdata = dict()
   xdata["frames"] = nf
   for k,v in zip(labels, aveu[1:]):
     xdata[k] = v
   cluster_energies[ic] = xdata

print cluster_energies
cpk.dump(cluster_energies,open(output,'w'),True)

