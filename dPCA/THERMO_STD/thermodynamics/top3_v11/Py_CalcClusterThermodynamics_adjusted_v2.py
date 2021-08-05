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


#thermodynamics_fe["-TdS"] = thermodynamics_fe["FE"] - thermodynamics_enthalpy["PE"]
#labels_fe.append("-TdS")

#thermodynamics_fe["-TdSsolv"] = thermodynamics_fe["-TdS"] - thermodynamics_fe["-TdSconf"]
#labels_fe.append("-TdSsolv")

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

def WriteThermodynamics_enthalpy_silent(thermodynamics, label):
  m = thermodynamics[label][0]
  sem = thermodynamics[label][1]
  return m, sem

def WriteThermodynamics_fe_silent(thermodynamics, label):
  m = np.mean(thermodynamics[label])
  sem = stats.sem(thermodynamics[label])
  return m, sem

for label in labels_enthalpy: WriteThermodynamics_enthalpy(thermodynamics_enthalpy, label)
for label in labels_fe: WriteThermodynamics_fe(thermodynamics_fe, label)

print "%10s"%"-TdS",
m_FE, sem_FE = WriteThermodynamics_fe_silent(thermodynamics_fe, "FE")
m_PE, sem_PE = WriteThermodynamics_enthalpy_silent(thermodynamics_enthalpy, "PE")
m_S = m_FE - m_PE
sem_S = np.sqrt(sem_FE**2+sem_PE**2) 
print "%8.3f%8.3f"%(m_S,sem_S)

print "%10s"%"-TdSsolv",
m_SC, sem_SC = WriteThermodynamics_fe_silent(thermodynamics_fe, "-TdSconf")
m_SS = m_S - m_SC
sem_SS = np.sqrt(sem_S**2+sem_SC**2)
print "%8.3f%8.3f"%(m_SS,sem_SS)


