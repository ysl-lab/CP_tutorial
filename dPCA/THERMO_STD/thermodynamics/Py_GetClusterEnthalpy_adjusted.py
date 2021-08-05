#!/bin/env python
import numpy as np



def flatten(l): return flatten(l[0]) + (flatten(l[1:]) if len(l) > 1 else []) if type(l) is list else [l]

def LoadEnthalpyInfoFile(inp):
  lines = [line.split() for line in open(inp,'r').readlines() if not line.startswith("#")]
  return flatten(lines)


def LoadEnthalpyFile (inp):
  lines = [line.split() for line in open(inp,'r').readlines() if not line.startswith("#")]
  lines = zip(*lines)
  return lines


import sys
inp1 = sys.argv[1]
inp2 = sys.argv[2]
icluster = int(sys.argv[3])+1

enthalpyComponents = LoadEnthalpyInfoFile("Enthalpy.info") 
enthalpyValues = LoadEnthalpyFile(inp1)
enthalpyValuesStd = LoadEnthalpyFile(inp2)

H1 = np.array(enthalpyValues[icluster], dtype=float)
H0 = np.array(enthalpyValues[-1], dtype=float)
dH = H1 - H0
H1_std = np.array(enthalpyValuesStd[icluster], dtype=float)
H0_std = np.array(enthalpyValuesStd[-1], dtype=float)
dH_std=np.sqrt(H1_std**2 + H0_std**2)
for ic in enthalpyComponents:
  idx1 = enthalpyValues[0].index(ic)
  print "%12s"%ic,
  print "%10.3f"%dH[idx1],
  print "%10.3f"%dH_std[idx1]
  

