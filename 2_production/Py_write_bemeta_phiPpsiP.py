#!/bin/bash

"""
python Py_write_bemeta_phiPpsiP.py prot.gro

Created by Diana Slough
Created on 06-26-2016

This program writes out 2D-BEMETA files that use phi' [H N CA C] and psi' [N CA C O]
NOTE: assumes CP and currently supports ALA and NMA...
"""

import sys
import numpy as np

NNEUTRAL=5    # number of neutral replicas

###############################################################

def read_lines(f):
   lines=[]
   for l in f: lines.append(l)
   return lines

def read_gro(lines):
   num_res=0
   resids=[]
   residues=[]
   atoms=[]
   indexes=[]
   for l in range(len(lines)):
      if l != 0 and l != 1:
         if l < (len(lines)-1):
            if int(lines[l].strip().split()[0][0]) > num_res: num_res=int(lines[l].strip().split()[0][0])
            resids.append(int(lines[l].strip().split()[0][0])) 
            residues.append(str(lines[l].strip().split()[0][1:]))
            atoms.append(str(lines[l].strip().split()[1]))
            indexes.append(int(lines[l].strip().split()[2]))
         else:  # box line
            break
   #print num_res,resids,residues,atoms,indexes
   return num_res,resids,residues,atoms,indexes

def get_phiP(num_res,resids,residues,atoms,indexes):
   phiP=[]
   for i in range(1,num_res+1):
      temp=[]
      for n in range(len(residues)):
         if resids[n] == i and residues[n] == 'ALA' and atoms[n] == 'H': temp.append(indexes[n])
         elif resids[n] == i and residues[n] == 'NMA' and atoms[n] == 'CN': temp.append(indexes[n])
      for n in range(len(residues)):
         if resids[n] == i and atoms[n] == 'N': temp.append(indexes[n])
      for n in range(len(residues)):
         if resids[n] == i and atoms[n] == 'CA': temp.append(indexes[n])
      for n in range(len(residues)):
         if resids[n] == i and atoms[n] == 'C': temp.append(indexes[n])
      phiP.append(temp)
   return phiP

def get_psiP(num_res,resids,residues,atoms,indexes):
   psiP=[]
   for i in range(1,num_res+1):
      temp=[]
      for n in range(len(residues)):
         if resids[n] == i and atoms[n] == 'N': temp.append(indexes[n])
      for n in range(len(residues)):
         if resids[n] == i and atoms[n] == 'CA': temp.append(indexes[n])
      for n in range(len(residues)):
         if resids[n] == i and atoms[n] == 'C': temp.append(indexes[n])
      for n in range(len(residues)):
         if resids[n] == i and atoms[n] == 'O': temp.append(indexes[n])
      psiP.append(temp)
   return psiP

###############################################################

def write_bemeta_files(num_res,phiP,psiP):
   num_replicas=num_res*2+NNEUTRAL
   for i in range(num_replicas):
      with open('bemeta.dat.%d' % (i),'w') as f:
         f.write('RANDOM_EXCHANGES\n\n')
         counter=0
         rep_counter=0
         for j in range(num_res):
            f.write('#Rep%d; Res%d:phiP/psiP\n' %(rep_counter,j+1))
            f.write('cv%d: TORSION ATOMS=%d,%d,%d,%d\n' % (counter,phiP[j][0],phiP[j][1],phiP[j][2],phiP[j][3]))
            f.write('cv%d: TORSION ATOMS=%d,%d,%d,%d\n\n' % (counter+1,psiP[j][0],psiP[j][1],psiP[j][2],psiP[j][3]))
            counter+=2
            rep_counter+=1 
         for j in range(num_res):
            if j != num_res-1:
               f.write('#Rep%d; Res%d:phiP/Res%d:psiP\n' %(rep_counter,j+1,j+2))
               f.write('cv%d: TORSION ATOMS=%d,%d,%d,%d\n' % (counter,psiP[j][0],psiP[j][1],psiP[j][2],psiP[j][3]))
               f.write('cv%d: TORSION ATOMS=%d,%d,%d,%d\n\n' % (counter+1,phiP[j+1][0],phiP[j+1][1],phiP[j+1][2],phiP[j+1][3]))
            else:
               f.write('#Rep%d; Res%d:phiP/Res%d:psiP\n' %(rep_counter,j+1,1))
               f.write('cv%d: TORSION ATOMS=%d,%d,%d,%d\n' % (counter,psiP[j][0],psiP[j][1],psiP[j][2],psiP[j][3]))
               f.write('cv%d: TORSION ATOMS=%d,%d,%d,%d\n\n' % (counter+1,phiP[0][0],phiP[0][1],phiP[0][2],phiP[0][3]))
            counter+=2
            rep_counter+=1 
         cvA=2*i
         cvB=2*i+1
         if i < num_res*2:
            f.write('METAD ARG=cv%d,cv%d SIGMA=0.31416,0.31416 HEIGHT=0.1 PACE=2000 LABEL=metad.cv%d FILE=HILLS\n' % (cvA,cvB,cvA) )
            f.write('PRINT ARG=cv%d,cv%d STRIDE=500 FILE=COLVAR\n\n' % (cvA,cvB) )
         else:
            f.write('#METAD ARG=cv%d,cv%d SIGMA=0.31416,0.31416 HEIGHT=0.1 PACE=2000 LABEL=metad.cv%d FILE=HILLS\n' % (cvA,cvB,cvA) )
            f.write('PRINT ARG=cv0 STRIDE=500 FILE=COLVAR\n\n') 
         f.write('ENDPLUMED\n')
   return

###############################################################

def main():
   with open(sys.argv[1],'r') as f:
      lines=read_lines(f)
   
   num_res,resids,residues,atoms,indexes=read_gro(lines)
   
   phiP=get_phiP(num_res,resids,residues,atoms,indexes)
   psiP=get_psiP(num_res,resids,residues,atoms,indexes)
   #print 'Phi',phiP
   #print 'Psi',psiP

   write_bemeta_files(num_res,phiP,psiP)

if __name__ == '__main__':
   main()