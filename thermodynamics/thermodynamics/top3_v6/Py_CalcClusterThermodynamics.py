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


labels = []
thermodynamics = dict()

for line in enthalpyValues:
  lab = line[0]
  labels.append(lab)
  val = np.array(line[1:], dtype=float)
  thermodynamics[lab] = val

for line in freeEnergyValues:
  lab = line[0]
  labels.append(lab)
  val = np.array(line[1:], dtype=float)
  thermodynamics[lab] = val


thermodynamics["-TdS"] = thermodynamics["FE"] - thermodynamics["PE"]
labels.append("-TdS")

thermodynamics["-TdSsolv"] = thermodynamics["-TdS"] - thermodynamics["-TdSconf"]
labels.append("-TdSsolv")

thermodynamics["EE_P"] = thermodynamics["EE14_P"] + thermodynamics["EEsr_P"] + thermodynamics["EElr_P"]
labels.append("EE_P")

thermodynamics["LJ_P"] = thermodynamics["LJ14_P"] + thermodynamics["LJsr_P"] + thermodynamics["LJlr_P"]
labels.append("LJ_P")

thermodynamics["EE_W"] = thermodynamics["EEsr_W"] + thermodynamics["EElr_W"]
labels.append("EE_W")

thermodynamics["LJ_W"] = thermodynamics["LJsr_W"] + thermodynamics["LJlr_W"]
labels.append("LJ_W")

thermodynamics["EE_PW"] = thermodynamics["EEsr_PW"] + thermodynamics["EElr_PW"]
labels.append("EE_PW")

thermodynamics["LJ_PW"] = thermodynamics["LJsr_PW"] + thermodynamics["LJlr_PW"]
labels.append("LJ_PW")

########################### OUTPUT ###########################

def WriteThermodynamics(thermodynamics, label):
  print "%10s"%label, 
  for x in thermodynamics[label]:
    print "%8.3f"%x,
  m = np.mean(thermodynamics[label])
  sem = stats.sem(thermodynamics[label])
  print "%8.3f%8.3f"%(m,sem)


for label in labels: WriteThermodynamics(thermodynamics, label)
