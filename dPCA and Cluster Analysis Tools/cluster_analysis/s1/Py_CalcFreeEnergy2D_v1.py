#!/usr/env python

import numpy as np
import sys

inp = sys.argv[1]
NX = int(sys.argv[2])

TEMPERATURE=300    # Kelvin

def CalcFreeEnergy (densit, temperature):
  kb = 0.008314           # Boltzmann constant kJ/mol/K
  kbt = kb*temperature
  fes = []
  for d in densit:
    if d <= 0.0: fes.append(10000.0)
    else: fes.append(-kbt*np.log(d))

# Shift the minimal value to 0
  fes = np.array(fes, dtype=float)
  fes -= np.min(fes) 
  return fes


lines = np.loadtxt(inp, dtype=np.float64, comments='#', usecols=(0,1,2), unpack=False)
densit = lines[:,2]
fes = CalcFreeEnergy(densit, TEMPERATURE)

lines[:,2] = fes
i = 1
for line in lines:
  x, y, z = line
  print x, y, z
  if i % NX == 0: print
  i += 1
  
