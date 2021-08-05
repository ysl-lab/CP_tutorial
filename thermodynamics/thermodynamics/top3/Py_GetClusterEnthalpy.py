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
inp = sys.argv[1]
ref = sys.argv[2]

enthalpyComponents = LoadEnthalpyInfoFile("Enthalpy.info") 
enthalpyValues = LoadEnthalpyFile(inp)
enthalpyValuesRef = LoadEnthalpyFile(ref)

H1 = np.array(enthalpyValues[1:], dtype=float)
H0 = np.array(enthalpyValuesRef[1:], dtype=float)
dH = H1 - H0
for ic in enthalpyComponents:
  idx1 = enthalpyValues[0].index(ic)
  idx2 = enthalpyValuesRef[0].index(ic)
  assert idx1 == idx2
  print "%12s"%ic,
  for x in dH[:,idx1]:
    print "%10.3f"%x,
  print 
  

