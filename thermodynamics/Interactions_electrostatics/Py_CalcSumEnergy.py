#1/bin/env python

import numpy as np
import sys
from scipy import stats

replicaID = int(sys.argv[1])  # replica ID
clusterID = int(sys.argv[2])


data = np.zeros((12,12)) 
interactionTypes = ["LJ-14", "LJ-SR", "Coul-14", "Coul-SR"]
for it in interactionTypes: 
  inp = "MAT_rep%d_cluster%d_%s.txt"%(replicaID, clusterID, it)
  lines = np.loadtxt(inp, comments="#", dtype=float)
  data += lines

out = "MAT_ENERGY_ALL_rep%d_cluster%d.txt"%(replicaID, clusterID) 
ofp = open(out, 'w')
for m1 in data:
  for m2 in m1:
    ofp.write("%8.3f "%m2)
  ofp.write("\n")
ofp.close()
