#!/usr/bin/python

"""
dPCA_project_onto_pc123.py phi_psi.xvg out_eigenvec.out
created by: Julia Rogers
created on: 10-12-15
modified on: 2-1-16

"""

import sys, os
import numpy as np
from numpy import matrix
from general import get_files_from_dir, get_lines
from math import cos, sin, radians
from CP_turns import get_only_phi_psi

def get_eigvecs():
   #TODO
   lines=get_lines(sys.argv[2])
   num=int(lines[1].split()[1])
   eigvecs=[]
   for i in range(7, len(lines), 7+num):
      temp=lines[i].split('{')[1].strip().strip('}')
      eigvecs.append([float(num.strip()) for num in temp.split(',')])
      for j in range(i+1, i+8):
         temp=lines[j].split('{')[1].strip().strip('}')
         eigvecs[-1].extend([float(num.strip()) for num in temp.split(',')])

   return eigvecs[1], eigvecs[2], eigvecs[3], eigvecs[4]

def calc_varience(frames, avs):
   var=[[frames[i][j]-avs[j] for j in range(len(avs))] for i in range(len(frames))]
   #print frames[0]
   #print avs
   #print var[0]
   return var

def project(frames, ev1, ev2, ev3, fname):
   with open(fname, 'w') as f:
      RT=matrix([ev1, ev2, ev3])
      for frame in frames:
         temp=matrix(frame)
         x=np.transpose(temp)
         #print x
         pcs=RT*x
         #print pcs
         #sys.exit()
         f.write('%15.5f%15.5f%15.5f\n' % (pcs[0], pcs[1], pcs[2]))

def main():
   av, ev1, ev2, ev3=get_eigvecs()
   #print av

   name=os.path.basename(sys.argv[1]).split('.')[0] #+'_'+os.path.basename(sys.argv[1]).split('_')[1]+'_'+os.path.basename(sys.argv[1]).split('_')[2]
   diheds=get_only_phi_psi(sys.argv[1])
   cossin=[]
   for dihed in diheds:
      curr=[]
      for d in dihed:
         curr.append(cos(radians(d)))
         curr.append(sin(radians(d)))
      for d in dihed:
         curr.append(cos(radians(d)))
         curr.append(sin(radians(d)))
      cossin.append(curr)

   var=calc_varience(cossin, av)
   project(var, ev1, ev2, ev3, name+'_pc1pc2pc3.txt')

if __name__ == '__main__':
   main()


