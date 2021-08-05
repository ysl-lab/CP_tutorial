#!/bin/env python

import numpy as np
import sys

inp = sys.argv[1]
out = sys.argv[2]
cid = int(sys.argv[3])


frameIDs   = np.loadtxt(inp,dtype=int,usecols=(0,),unpack=False)
clusterIDs = np.loadtxt(inp,dtype=int,usecols=(-1,),unpack=False)
clusterFrames = frameIDs[np.where(clusterIDs == cid)]

nframes = len(frameIDs)
ncframes = len(clusterFrames)
Pop = float(ncframes)/float(nframes)
ofp = open(out,'w')
ofp.write("[ state %d, Pop %8.3f ]\n"%(cid, Pop*100))
for ifm in range(ncframes):
  ofp.write("%10d"%clusterFrames[ifm])
  if (ifm+1)%10 == 0: ofp.write("\n") 

ofp.write("\n")

print "%8d %8.3f"%(cid, Pop*100)
