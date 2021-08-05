#1/bin/env python

import numpy as np
import sys
from scipy import stats

clusterID = int(sys.argv[1])


data = []


for it in range(12,17):
  inp = "MAT_ENERGY_ALL_rep%d_cluster%d.txt"%(it, clusterID)
  lines = np.loadtxt(inp, comments="#", dtype=float)
  data.append(lines)
 
data = np.array(data, dtype=float)
mean = np.mean(data, axis=0)
sem = stats.sem(data, axis=0, ddof=1)

out = "MAT_MEAN_cluster%d.txt"%(clusterID)
ofp = open(out, 'w')
for m1 in mean:
  for m2 in m1:
    ofp.write("%8.3f "%m2)
  ofp.write("\n")
ofp.close()

out = "MAT_SEM_cluster%d.txt"%(clusterID)
ofp = open(out, 'w')
for m1 in sem:
  for m2 in m1:
    ofp.write("%8.3f "%m2)
  ofp.write("\n")
ofp.close()
  

