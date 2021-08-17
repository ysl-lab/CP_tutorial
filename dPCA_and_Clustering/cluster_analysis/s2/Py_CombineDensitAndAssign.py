#!/bin/env python

import numpy as np
import sys

inp1 = sys.argv[1]
inp2 = sys.argv[2]

lines1 = np.loadtxt(inp1, dtype=np.float, usecols=(0,1,2,3),unpack=False)
lines2 = np.loadtxt(inp2, dtype=np.int, usecols=(1,),unpack=False)
n = len(lines1) 
for i in range(n):
  x, y, z, d = lines1[i]
  a = lines2[i]
  if a != 0:
    print "%10.5f %10.5f %10.5f %5d %10.5f"%(x, y, z, a, d)
