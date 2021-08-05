#!/bin/env python
#-*- coding: utf-8 -*-


import numpy as np
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
from Py_GenTableGlobalVar  import *

import sys
tbl=sys.argv[1]


def LoadData (inp):
  labels = [line.split('"')[1] for line in open(inp, 'r').readlines() if not line.startswith('#')]
  data = [line.split('"')[2] for line in open(inp, 'r').readlines() if not line.startswith('#')]
  mean = [i.split()[-2] for i in data]
  sem  = [i.split()[-1] for i in data]
  return labels, mean , sem


allMeans = []
allSEMs = []
for ic in reversed(range(IC,NC+1)):
  label,mean, sem  = LoadData('%s_cluster%d.dat'%(tbl,ic))
  print label
  allMeans.append(mean)
  allSEMs.append(sem)

allMeans = np.array(allMeans, dtype=float)
allSEMs = np.array(allSEMs, dtype=float)
print allMeans

NS = len(allMeans)
xticks = label
yticks = ["state %d"%i for i in range(1,NS+1)]
print xticks

fig = plt.figure(figsize=(9.8,2.0))
rect = [0.02,0.04,0.95,0.75] # l, b, w ,h
ax = fig.add_axes(rect)
cax=ax.matshow(allMeans, aspect='auto', cmap='bwr', vmin=vmin, vmax=vmax, interpolation='none')
cbar = fig.colorbar(cax,ticks=range(vmin,vmax+1,vstep), pad=0.01)  # pad is the distance between colorbar and the table
ax.annotate("kJ/mol", xy=(0.88,0.5), xycoords='figure fraction', ha='center', va='center', rotation='vertical',size=12)

ax.tick_params(labelsize=16)
#ax.set(xticks=range(len(xticks)), xticklabels=xticks, yticks=range(len(yticks)), yticklabels=yticks)
ax.set(xticks=range(len(xticks)), xticklabels=xticks, yticks=range(len(yticks)), yticklabels=[])

# Write the numbers 
for (i,j), val in np.ndenumerate(allMeans):
  color="black"
  if val < 0.75*vmin or val > 0.75*vmax: color="white"
  if i == 0:
    val = 0.0
    ax.annotate("%6.2f"%(val,), (j,i), ha='center', va='center', size=12, color=color,weight='bold')
  else:
    ax.annotate("%6.2f$\pm$%4.2f"%(val,allSEMs[i,j]), (j,i), ha='center', va='center', size=12, color=color,weight='bold')


fig.savefig("%s_Thermodynamics.png"%tbl)

