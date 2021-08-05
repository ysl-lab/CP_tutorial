#!/bin/python 

"""
python Py_average_populations.py .

Created by Diana Slough
Created on 07-11-2016

Averages populations and std is standard error of means
"""

import sys,os
import fnmatch
import numpy as np
from scipy import stats

def read_lines(f1):
   """
   Takes open file 
   Returns list of each line in the file
   """
   lines=[]
   for l in f1:
      lines.append(l)
   return lines

def write_output(clusters):
   with open('Out_average_populations.txt','w') as f:
      f.write('# cluster average sem std\n')
      for k in sorted(clusters.keys()):
         ave=np.mean(clusters[k])
         sem=stats.sem(clusters[k])
         std=np.std(clusters[k])
         f.write('%d %8.3f %8.3f %8.3f\n' % (k,ave,sem,std))
   return

def main():
   files=[]
   for file in os.listdir(sys.argv[1]):   # looks for file names in the directory specified by first command line arguement && appends to list of files
      if fnmatch.fnmatch(file, 'population*.txt'): files.append(file)

   clusters=dict()
   for f in files:
      with open(f,'r') as f1:
         lines=read_lines(f1)
      
         for l in lines:
            if int(l.strip().split()[0]) not in clusters.keys():
               clusters[int(l.strip().split()[0])]=[]
               clusters[int(l.strip().split()[0])].append(float(l.strip().split()[-1]))
            else: clusters[int(l.strip().split()[0])].append(float(l.strip().split()[-1]))
   print clusters

   write_output(clusters)

if __name__ == '__main__': 
   main()
