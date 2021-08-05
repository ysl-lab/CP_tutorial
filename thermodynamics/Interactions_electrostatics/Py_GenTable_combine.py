#!/bin/env python
#-*- coding: utf-8 -*-


import numpy as np
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colors

import sys
inp_m1=sys.argv[1]
inp_e1=sys.argv[2]
inp_m2=sys.argv[3]
inp_e2=sys.argv[4]
out = sys.argv[5]

def genColorMap(cmap):
# Creat new colormap with white for 1.0
#  cvals = [(cmap(i)) for i in xrange(0,255)] + [('white')]
  cvals = [('gray')] + [(cmap(i)) for i in xrange(1,255)] + [("gray")]
#  cvals = [(cmap(i)) for i in xrange(1,256)] 
  new_map = colors.LinearSegmentedColormap.from_list('new_map',cvals, N=256)
  return new_map


mean1  = np.loadtxt(inp_m1, dtype=float, comments="#")
sem1   = np.loadtxt(inp_e1, dtype=float, comments="#")
mean2  = np.loadtxt(inp_m2, dtype=float, comments="#")
sem2   = np.loadtxt(inp_e2, dtype=float, comments="#")

label = ["NH_1", "CO_1", "NH_2", "CO_2", "NH_3", "CO_3", "NH_4", "CO_4", "NH_5", "CO_5", "NH_6", "CO_6",] 

mean = np.empty_like(mean1)
mean[:] = mean2
sem = np.empty_like(sem1)
sem[:] = sem2
n = len(label)
for i in range(n):
  for j in range(i+1,n):
    mean[i,j] = mean1[i][j]
    sem[i,j] = sem1[i][j]
  


vmin=np.min(mean)
vmax=np.max(mean)
#vmin=-80
#vmax=80
if np.abs(vmin) > np.abs(vmax):
  vmax = np.abs(vmin)
else:
  vmin = -np.abs(vmax)

NS = len(mean)
xticks = label
yticks = label

fig = plt.figure(figsize=(10,4.0))
rect = [0.09,0.08,0.999,0.8] # l, b, w ,h
ax = fig.add_axes(rect)
#ax = fig.add_subplot(111)
#ax.matshow(allMeans, aspect='auto', cmap='bwr')
cmap=genColorMap(cmx.bwr)
#cax=ax.matshow(mean, aspect='auto', cmap='bwr', vmin=vmin, vmax=vmax, interpolation='none')
cax=ax.matshow(mean, aspect='auto', cmap=cmap, vmin=-20, vmax=20, interpolation='none')
#cbar = fig.colorbar(cax,ticks=range(-80,81,40), pad=0.01)  # pad is the distance between colorbar and the table
cbar = fig.colorbar(cax,ticks=range(-20,21,10), pad=0.01)  # pad is the distance between colorbar and the table
ax.annotate("kJ/mol", xy=(0.98,0.5), xycoords='figure fraction', ha='center', va='center', rotation='vertical',size=12)
ax.autoscale(False)
ax.plot(range(-1,NS+2), range(-1,NS+2), lw=4, c="black")

ax.tick_params(labelsize=12)
ax.set(xticks=range(len(xticks)), xticklabels=xticks, yticks=range(len(yticks)), yticklabels=yticks)

# Write the numbers 
for (i,j), val in np.ndenumerate(mean):
   ax.annotate("%5.1f$\pm$%2.1f"%(val,sem[i,j]), (j,i), ha='center', va='center', size=9)

fig.savefig(out)

