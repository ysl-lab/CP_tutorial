#!/bin/python

"""
python Py_write_dPCA_assign_fortran.py all.txt 

Created by Diana Slough
Created on 1-5-2017

This program reads all.txt (from dPCA) and writes Assign.f90
"""

import sys,os
import numpy as np
import math

##########################################################

def read_lines(f):
   lines=[]
   for l in f: lines.append(l)
   return lines

def get_min_max(lines):
   x=[]
   y=[]
   z=[]
   for l in lines:
      spl=l.strip().split()
      x.append(float(spl[0]))
      y.append(float(spl[1]))
      z.append(float(spl[2]))
   min_max=[]
   min_max.append(math.floor(min(x)))
   min_max.append(math.ceil(max(x)))
   min_max.append(math.floor(min(y)))
   min_max.append(math.ceil(max(y)))
   min_max.append(math.floor(min(z)))
   min_max.append(math.ceil(max(z)))
   #print min_max
   return min_max 

##########################################################

def write_output(min_max):
   with open('Assign.f90','w') as f:
      f.write('program assign\n')
      f.write('  implicit none\n')
      f.write('  integer, parameter :: ng = 50*50*50\n')
      f.write('  integer, dimension (ng) :: asn\n')
      f.write('  double precision, dimension(ng):: gx, gy, gz\n')
      f.write('  double precision:: x, y, z, xsize_half, ysize_half, zsize_half\n')
      f.write('  integer :: ng_kept,ig, ig_max, cl, ik\n')
      f.write('  integer :: ios\n')
      f.write('  character(len=1000):: input, output\n\n')
      f.write('  gx = 0.d0\n')
      f.write('  gy = 0.d0\n')
      f.write('  gz = 0.d0\n')
      f.write('  asn = 0.d0\n\n')
      f.write('  call get_command_argument(1,input)\n')
      f.write('  call get_command_argument(2,output)\n\n')
      f.write('  xsize_half = 0.5*(%.1f+%.1f)/50.0\n' % (np.abs(min_max[0]),min_max[1]))
      f.write('  ysize_half = 0.5*(%.1f+%.1f)/50.0\n' % (np.abs(min_max[2]),min_max[3]))
      f.write('  zsize_half = 0.5*(%.1f+%.1f)/50.0\n\n' % (np.abs(min_max[4]),min_max[5]))
      f.write('! Read the assignments for the density grid\n')
      f.write('  open(12,file=\'../GRID_ASSIGNATION2\',action=\'read\')\n')
      f.write('  ig = 1\n')
      f.write('  do\n')
      f.write('    read(12,*, iostat=ios) gx(ig), gy(ig), gz(ig), asn(ig)\n')
      f.write('    if (ios /= 0) exit\n')
      f.write('    ig = ig + 1\n')
      f.write('  end do\n')
      f.write('  close(12)\n')
      f.write('  ig_max = ig\n\n')
      f.write('  write(6,*) \'Number of assigned grids: \', ig_max\n\n')
      f.write('  open(12,file=input, action=\'read\')\n')
      f.write('  open(13,file=output, action=\'write\')\n\n')
      f.write('  ik = 1\n')
      f.write('  do\n')
      f.write('    read(12,*,iostat=ios) x, y, z\n')
      f.write('    if (ios /= 0) exit\n')
      f.write('    cl = 0\n')
      f.write('    do ig = 1, ig_max\n')
      f.write('      if (abs(x-gx(ig)) < xsize_half .and. abs(y-gy(ig)) < ysize_half .and. abs(z-gz(ig)) < zsize_half) then\n')
      f.write('        cl = asn(ig)\n')
      f.write('      end if\n')
      f.write('    end do\n')
      f.write('    write(13,\'(i8,f10.5,f10.5,f10.5,i8)\') ik, x, y, z, cl\n')
      f.write('    ik = ik+1\n')
      f.write('  end do\n')
      f.write('  close(12)\n')
      f.write('  close(13)\n\n')
      f.write('end program assign\n')
   return

##########################################################

def main():
   with open(sys.argv[1],'r') as f:
      lines=read_lines(f)
   min_max=get_min_max(lines)
   write_output(min_max)

if __name__ == '__main__': 
   main()
