#!/bin/env python
#-*- coding: utf-8 -*-


import numpy as np
import matplotlib
matplotlib.use('AGG')
from Py_GenTableGlobalVar  import *
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
import sys
tbl=sys.argv[1]


def LoadData (inp):
  labels = [line.split('"')[1] for line in open(inp, 'rb').readlines() if not line.startswith('#')]
  data = [line.split('"')[2] for line in open(inp, 'r').readlines() if not line.startswith('#')]
  mean = [i.split()[-2] for i in data]
  sem  = [i.split()[-1] for i in data]
  return labels, mean , sem

def GenColor (inp):
  data = [line.split('"')[2] for line in open(inp, 'r').readlines() if not line.startswith('#')]
  mean_c = [i.split()[-2] if (float(i.split()[-2]) > 0 and (float(i.split()[-2]) - float(i.split()[-1])) > 0) or (float(i.split()[-2]) < 0 and (float(i.split()[-2]) + float(i.split()[-1])) < 0) else "0.0" for i in data]
  return mean_c


allMeans = []
allSEMs = []
meanColor = []
for ic in reversed(range(IC,NC+1)):
  label,mean, sem  = LoadData('%s_cluster%d.dat'%(tbl,ic))
  mean_color = GenColor('%s_cluster%d.dat'%(tbl,ic))
#  print mean_color
#  print label
  allMeans.append(mean)
  allSEMs.append(sem)
  meanColor.append(mean_color)

allMeans = np.array(allMeans, dtype=float)
allSEMs = np.array(allSEMs, dtype=float)
meanColor = np.array(meanColor, dtype=float)
#print allMeans

NS = len(allMeans)
xticks = label
yticks = ["state %d"%i for i in range(1,NS+1)]
#print xticks

fig = plt.figure(figsize=(4.9,2.0))
rect = [0.02,0.04,0.95,0.70] # l, b, w ,h
ax = fig.add_axes(rect)
cax=ax.matshow(meanColor, aspect='auto', cmap='bwr', vmin=vmin, vmax=vmax, interpolation='none')
cbar = fig.colorbar(cax,ticks=range(vmin,vmax+1,vstep), pad=0.03)  # pad is the distance between colorbar and the table
cbar.ax.tick_params(labelsize=18)
ax.annotate("\pmb{kJ/mol}", xy=(0.96,0.5), xycoords='figure fraction', ha='center', va='center', rotation='vertical',size=18, weight='bold')

ax.tick_params(labelsize=22)
#ax.set(xticks=range(len(xticks)), xticklabels=xticks, yticks=range(len(yticks)), yticklabels=yticks)
ax.set(xticks=range(len(xticks)), xticklabels=xticks, yticks=range(len(yticks)), yticklabels=[])

# Write the numbers 
for (i,j), val in np.ndenumerate(allMeans):
  color="black"
#  if val < 0.75*vmin or val > 0.75*vmax: color="white"
  if i == 0:
    val = 0.0
    ax.annotate("\pmb{%6.2f}"%(val,), (j,i), ha='center', va='center', size=16, color=color,weight='bold')
  else:
    ax.annotate("\pmb{%6.2f$\pm$%4.2f}"%(val,allSEMs[i,j]), (j,i), ha='center', va='center', size=16, color=color,weight='bold')

fig.savefig("%s_Thermodynamics.png"%tbl)

