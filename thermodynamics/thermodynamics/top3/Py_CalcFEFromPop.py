#!/bin/python

import numpy as np
import cPickle as cpk

#kb=0.0019872041    # kcal/mol/K
kb=0.0083144621     # kJ/mol/K
T=300              # K


def calcFreeEnergyFromPop(pop,maxP, T):
   global kb
   dFE=-kb*T*np.log(pop/maxP)
   return dFE 

import sys
popFile = sys.argv[1]
ic = int(sys.argv[2])

pops       = np.loadtxt(popFile, dtype=float, usecols=(-1,))
clusterIDs = np.loadtxt(popFile, dtype=int, usecols=(0,))

maxP = np.max(pops)
n = len(clusterIDs)

for i in range(n):
  if clusterIDs[i] ==  ic:
    x = calcFreeEnergyFromPop(pops[i], maxP, 300.0)
    print "%10.3f"%x
    

