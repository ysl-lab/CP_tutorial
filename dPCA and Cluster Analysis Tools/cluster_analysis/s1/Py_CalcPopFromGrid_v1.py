#!/bin/env python

import numpy as np
import sys


xmin = float(sys.argv[1])
xmax = float(sys.argv[2])
ymin = float(sys.argv[3])
ymax = float(sys.argv[4])
zmin = float(sys.argv[5])
zmax = float(sys.argv[6])
nx   = int(sys.argv[7])
ny   = int(sys.argv[8])
nz   = int(sys.argv[9])


area_tot = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)
xsize = (xmax - xmin)/nx
ysize = (ymax - ymin)/ny
zsize = (zmax - zmin)/nz
area_grid = xsize*ysize*zsize

states = np.loadtxt("GRID_ASSIGNATION", dtype=int,usecols=(3,))
densit = np.loadtxt("GRID_ASSIGNATION", dtype=float,usecols=(4,))

ns = max(states) 
nl = len(states)
print ">> Number of clusters: ", ns

pop = np.zeros(ns)
for i in range(nl):
  s = states[i] - 1
  pop[s] += densit[i]*area_grid

print "Populations calculated from grids (%):"
for i in range(ns):
  print "  cluster %d: %8.3f"%(i+1, pop[i]*100)

#################################################################
ordered_pop = sorted(range(len(pop)), key=lambda k: pop[k])
new_states = [ordered_pop.index(i)+1 for i in range(ns)]

lines1 = np.loadtxt('GRID_ASSIGNATION', usecols=(0,1,2), dtype=float)
lines2 = np.loadtxt('GRID_ASSIGNATION', usecols=(3,), dtype=int)

ofp = open('GRID_ASSIGNATION2', 'w')
n = len(lines1)
for i in range(n):
  x, y, z = lines1[i]
  s = lines2[i] - 1
  snew = new_states[s] 
  ofp.write("%8.3f %8.3f %8.3f %8d\n"%(x, y, z, snew))
ofp.close()
