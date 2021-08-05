#!/bin/env python
#-*- coding: utf-8 -*-


import numpy as np
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import sys
tbl=sys.argv[1]


IC = int(sys.argv[2])
NC = int(sys.argv[3])

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

vmin=np.min(allMeans)
vmax=np.max(allMeans)
vmin=-80
vmax=80
if np.abs(vmin) > np.abs(vmax):
  vmax = np.abs(vmin)
else:
  vmin = -np.abs(vmax)

NS = len(allMeans)
xticks = label
yticks = ["state %d"%i for i in range(1,NS+1)]
print xticks

fig = plt.figure(figsize=(10,3.0))
rect = [0.09,0.08,0.999,0.8] # l, b, w ,h
ax = fig.add_axes(rect)
#ax = fig.add_subplot(111)
#ax.matshow(allMeans, aspect='auto', cmap='bwr')
cax=ax.matshow(allMeans, aspect='auto', cmap='bwr', vmin=vmin, vmax=vmax, interpolation='none')
cbar = fig.colorbar(cax,ticks=range(-80,81,40), pad=0.01)  # pad is the distance between colorbar and the table
ax.annotate("kJ/mol", xy=(0.98,0.5), xycoords='figure fraction', ha='center', va='center', rotation='vertical',size=12)

ax.tick_params(labelsize=12)
ax.set(xticks=range(len(xticks)), xticklabels=xticks, yticks=range(len(yticks)), yticklabels=yticks)

# Write the numbers 
for (i,j), val in np.ndenumerate(allMeans):
  if i == 0:
    ax.annotate("%6.2f"%(val,), (j,i), ha='center', va='center', size=9)
  else:
    ax.annotate("%6.2f$\pm$%4.2f"%(val,allSEMs[i,j]), (j,i), ha='center', va='center', size=9)

fig.savefig("%s_Thermodynamics.png"%tbl)

