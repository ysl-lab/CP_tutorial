#!/bin/env


import sys

inp = sys.argv[1]

imps = [line.strip() for line in open("Impropers.dat", 'r')]
lines = [line for line in open(inp, 'r')]

flag = False
n = len(lines)
for i in range(n):
  if lines[i] == "; Include Position restraint file\n":
    a1, a2, a3, a4 = imps[0].split()
    xx = "%5s %5s %5s %5s     4     ; added by DPS"%(a1,a2,a3,a4)
    print xx
    a1, a2, a3, a4 = imps[1].split()
    xx = "%5s %5s %5s %5s     4     ; added by DPS"%(a1,a2,a3,a4)
    print xx
    flag = True
  if i > 0:
    print lines[i-1],

print lines[n-1],

if not flag:
  print "********************************************************"
  print "**ERROR: NO dihedral is added "
