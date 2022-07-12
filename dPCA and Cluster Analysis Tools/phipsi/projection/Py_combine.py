#!/bin/env python

import numpy as np
import sys

inp1 = sys.argv[1]
inp2 = sys.argv[2]
#Combines the 2D PC projections to create a 3D PC projection
lines1 = np.loadtxt(inp1, dtype=float, usecols=(0,1), unpack=False)
lines2 = np.loadtxt(inp2, dtype=float, usecols=(1,), unpack=False)

n = len(lines1)
for i in range(n):
   pc1,pc2 = lines1[i]
   pc3 = lines2[i]
   print pc1,pc2,pc3
