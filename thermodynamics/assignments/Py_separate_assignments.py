#!/bin/python 

"""
python Py_separate_assignments.py assignments.txt Out_num_indexes_corr_struct.txt traj_order.txt

Created by Diana Slough
Created on 8-20-2016

This program separates the assignments files into individual replicas
"""

import sys,os
import numpy as np

def read_lines(f):
   lines=[]
   for l in f: lines.append(l)
   return lines

def write_output(trj_order,new_assign):
   counter=0
   for i in trj_order:
      with open('assignment_rep%d.txt' % i,'w') as f:
         c2=1
         for i in new_assign[counter]:
            spl=i.strip().split()
            f.write('%8d %9.5f %9.5f %9.5f %7d\n' % (c2,float(spl[1]),float(spl[2]),float(spl[3]),int(spl[4]))) 
            c2+=1
      counter+=1
   return

def main():
   trj_order=[]
   with open(sys.argv[3],'r') as f:
      lines=read_lines(f)
      for l in lines:
         spl=l.strip().split()
         for i in spl: trj_order.append(int(i))
   print 'trajectory order:',trj_order
      
   num_fr=dict()
   with open(sys.argv[2],'r') as f:
      lines=read_lines(f)
      for l in lines:
         spl=l.strip().split()
         num_fr[int(str(spl[0][4])+str(spl[0][5]))]=int(float(spl[1]))
   print 'Number of frames:',num_fr

   new_assign=[]
   with open(sys.argv[1],'r') as f:
      lines=read_lines(f)
      counter=0
      for i in trj_order:
         temp=[]
         for j in range(num_fr[i]): 
            temp.append(lines[counter])
            counter+=1 
         new_assign.append(temp)
   write_output(trj_order,new_assign)

if __name__ == '__main__': 
   main()
