#!/bin/env python

import numpy as np
import sys

inp = sys.argv[1]

lines = np.loadtxt(inp,dtype=int,usecols=(0,4),unpack=False)
nstates = np.max(lines[:,1])+1

index = []
for s in range(nstates):
  index.append([])


for line in lines:
  iframe, st = line
  index[st].append(iframe)


i = 0
for idx in index:
  print "[ state_%d ]"%i
  ic = 0
  for ix in idx:
    print "%10d"%ix,
    ic += 1
    if ic%10 == 0: print
  print 
  i += 1

