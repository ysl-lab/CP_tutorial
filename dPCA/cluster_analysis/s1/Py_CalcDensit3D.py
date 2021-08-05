#!/usr/env python

import numpy as np
import sys

inp = sys.argv[1]
xmin = float(sys.argv[2])
xmax = float(sys.argv[3])
ymin = float(sys.argv[4])
ymax = float(sys.argv[5])
zmin = float(sys.argv[6])
zmax = float(sys.argv[7])
nx   = int(sys.argv[8])
ny   = int(sys.argv[9])
nz   = int(sys.argv[10])

area_tot = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)
xsize = (xmax - xmin)/nx
ysize = (ymax - ymin)/ny
zsize = (zmax - zmin)/nz
area_grid = xsize*ysize*zsize

# Load the (x, y) values
lines = np.loadtxt(inp, dtype=np.float64, comments='#', usecols=(0,1,2), unpack=False)

denst2d = np.zeros([nx, ny, nz])
for line in lines:
  x, y, z = line
  ix = int((x - xmin)/xsize)
  iy = int((y - ymin)/ysize)
  iz = int((z - zmin)/zsize)
  denst2d[ix,iy,iz] += 1

norm = np.sum(denst2d)
denst2d = denst2d/norm/area_grid

#print np.sum(denst2d*area_grid)

for ix in range(nx):
  for iy in range(ny):
    for iz in range(nz):
      x = xmin + (ix + 0.5) * xsize
      y = ymin + (iy + 0.5) * ysize
      z = zmin + (iz + 0.5) * zsize
      d = denst2d[ix, iy, iz]
      if d > 0.0:
         print "%-15.5f %-15.5f %-15.5f %-15.8f"%(x, y, z, d)
  
