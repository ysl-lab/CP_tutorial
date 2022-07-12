#!/bin/env python

import numpy as np
import math
import sys

inp = sys.argv[1]

lines = np.loadtxt(inp, dtype=np.float, comments='#', usecols=(0,1,2,3), unpack=False)
nlines = len(lines)

for i in range(nlines):
  for j in range(i+1, nlines):
    x0, y0, z0, d0 = lines[i]
    x1, y1, z1, d1 = lines[j]
#    dist = math.hypot(x1-x0, y1-y0, z1-z0)
    dist = math.sqrt((x1-x0)**2 + (y1-y0)**2 + (z1-z0)**2)
    print "%10d %10d %15.5f %10.5f %10.5f"%(i+1, j+1, dist, d0, d1)
