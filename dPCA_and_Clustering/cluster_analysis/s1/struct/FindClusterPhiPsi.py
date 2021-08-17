#!/bin/env python

import numpy as np
import sys

asn=sys.argv[1]
ang=sys.argv[2]
out=sys.argv[3]


def LoadDihedrals (inp): 
  lines = [i.split()[2:] for i in open(inp,'r').readlines() if not (i.startswith('#') or i.startswith('@'))]
  return np.array(lines,dtype=float)

dihedrals = LoadDihedrals(ang)
clusters  = np.loadtxt(asn,dtype=int,usecols=(3,),unpack=False)

nf1 = len(dihedrals)
nf2 = len(clusters)
if not nf1 == nf2: 
  print "**Error: Number of frames does not match ..."
  sys.exit(1)

ncl = np.max(clusters) + 1
print ">> Number of Clusters: ", ncl

dihedrals_clusters = []
for icl in range(ncl):
  dihedrals_clusters.append([])

for i in range(nf1):
  dihed = dihedrals[i]
  icl = clusters[i]
  dihedrals_clusters[icl].append(dihed)

  
for icl in range(ncl):
  dihedrals = dihedrals_clusters[icl]
  output="%s_cluster_%d.txt"%(out, icl)
  with open(output,'w') as ofp:
    for dihed in dihedrals:
      s  = " "
      for d in dihed:
        s += "%8.3f "%d
      s += "\n"
      ofp.write(s)

