#!/bin/env python


import numpy as np
import sys
import cPickle as cpk 
from scipy import stats

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


def CalcClusterEnergy (labels, energyInfo, energies):
  tol = 0.01
  #PE,PE2,Nonbonded,LJ,LJsr,EE,EEsr,EElr,PV
  hw=["PE","PE2","Nonbonded","LJ","LJsr","LJlr","EE","EEsr","EElr","PV"]
  clusterEnergies = dict()
  clusterEnergies_std = dict()
  forsave = dict()
  for key, vals in energyInfo.iteritems():
    ene = 0.0
    for v in vals: ene += energies[:,labels.index(v)+1]  #plus 1 because the first column is timestep
    if key in hw:
       forsave[key]=ene
    clusterEnergies[key] = np.mean(ene,axis=0)             # average energy
    clusterEnergies_std[key] = np.std(ene,axis=0,ddof=1)

  assert (clusterEnergies['PE'] - clusterEnergies['PE2']) < tol  # Check if PE and PE2 are idential within some tolence
  return clusterEnergies, clusterEnergies_std, forsave

#################################################################################################################
nrep=5
ipindex=2*nrep+1
energyConfigFile = sys.argv[ipindex]
output=sys.argv[ipindex+1]

u=dict()

for x in range(0,nrep):
   eindex=2*x+1
   aindex=2*x+2
   energyFile = sys.argv[eindex]
   assignmentFile = sys.argv[aindex]

#   print ">>Input Energy File: %s"%energyFile
#   print ">>Input Assignment File: %s"%assignmentFile

   labels   = GetLegend(energyFile)
   energies = GetEnergies(energyFile)
   cids = np.loadtxt(assignmentFile, dtype=int, usecols=(-1,))
   assert len(cids) == len(energies)

   nc = np.max(cids)  # number of clusters
   for ic in range(nc+1):
      if x == 0:
        u[ic] =  energies[np.where(cids==ic)]  # energy of each cluster
      else:
        u[ic] = np.concatenate((u[ic], energies[np.where(cids==ic)]), axis=0)

cluster_energies = dict()
cluster_energies_std = dict()
energyInfo, energyTerms = LoadEnergyConfigInfo(energyConfigFile)
print "#",
for k in energyTerms:
  print "%10s "%k,
print

save_energies = dict()
for ic in range(nc+1):
   aveu, stdu, savedata = CalcClusterEnergy(labels, energyInfo, u[ic])
   save_energies[ic]=savedata
   nf = len(u[ic])                        # number of frames in each cluster
#  print "Number of frames for cluster %d: %d"%(ic, nf)

   print "# ClusterID = %d average energy"%ic
   print " ",
   for k in energyTerms:
     print "%10.3f "%aveu[k],
   print

print "#",
for k in energyTerms:
  print "%10s "%k,
print

for ic in range(nc+1):
   aveu, stdu, temp = CalcClusterEnergy(labels, energyInfo, u[ic])
   nf = len(u[ic])                        # number of frames in each cluster
#  print "Number of frames for cluster %d: %d"%(ic, nf)

   print "# ClusterID = %d std"%ic
   print " ",
   for k in energyTerms:
     print "%10.3f "%stdu[k],
   print

cpk.dump(save_energies,open(output,'w'),True)
