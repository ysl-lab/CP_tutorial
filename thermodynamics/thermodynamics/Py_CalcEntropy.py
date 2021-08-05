#!/bin/python

import numpy as np
import cPickle as cpk

import sys
entFile = sys.argv[1]
ic = int(sys.argv[2])

entropies   = np.loadtxt(entFile, dtype=float, usecols=(-1,))
clusterIDs = np.loadtxt(entFile, dtype=int, usecols=(0,))

ent0 = entropies[-1]  # note: the data readed in is TdS, not -TdS
n = len(clusterIDs)
for i in range(n):
  if clusterIDs[i] ==  ic:
    xx = ent0 - entropies[i] 
    print "%10.3f"%xx
    

