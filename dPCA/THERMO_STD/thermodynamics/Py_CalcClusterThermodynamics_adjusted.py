#!/bin/env python
import numpy as np
from scipy import stats

def LoadEnthalpyFile (inp):
  lines = [line.split() for line in open(inp,'r').readlines() if not line.startswith("#")]
  return lines


import sys
enthalpyFile = sys.argv[1]
freeEnergyFile = sys.argv[2]

enthalpyValues = LoadEnthalpyFile(enthalpyFile)
freeEnergyValues = LoadEnthalpyFile(freeEnergyFile) 


labels_enthalpy = []
thermodynamics_enthalpy = dict()

labels_fe = []
thermodynamics_fe = dict()

for line in enthalpyValues:
  lab = line[0]
  labels_enthalpy.append(lab)
  val = np.array(line[1:], dtype=float)
  thermodynamics_enthalpy[lab] = val

for line in freeEnergyValues:
  lab = line[0]
  labels_fe.append(lab)
  val = np.array(line[1:], dtype=float)
  thermodynamics_fe[lab] = val


#thermodynamics["-TdS"] = thermodynamics["FE"] - thermodynamics["PE"]
#labels.append("-TdS")

#thermodynamics["-TdSsolv"] = thermodynamics["-TdS"] - thermodynamics["-TdSconf"]
#labels.append("-TdSsolv")

########################### OUTPUT ###########################

def WriteThermodynamics_enthalpy(thermodynamics, label):
  print "%10s"%label, 
  m = thermodynamics[label][0]
  sem = thermodynamics[label][1]
  print "%8.3f%8.3f"%(m,sem)

def WriteThermodynamics_fe(thermodynamics, label):
  print "%10s"%label,
  #for x in thermodynamics[label]:
  #  print "%8.3f"%x,
  m = np.mean(thermodynamics[label])
  sem = stats.sem(thermodynamics[label])
  print "%8.3f%8.3f"%(m,sem)

for label in labels_enthalpy: WriteThermodynamics_enthalpy(thermodynamics_enthalpy, label)
for label in labels_fe: WriteThermodynamics_fe(thermodynamics_fe, label)
