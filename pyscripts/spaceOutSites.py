#!/usr/bin/env python
import sys

def my_compare1st(x, y):
    if float(x[0]) > float(y[0]):
        return 1
    elif float(x[0])==float(y[0]):
        return 0
    else: # x<y
        return -1

def my_compare(x, y):
    if float(x[2]) < float(y[2]):
        return 1
    elif float(x[2])==float(y[2]):
        return 0
    else: # x<y
        return -1

def my_compare2(x, y):
    if float(x[3]) < float(y[3]):
        return 1
    elif float(x[3])==float(y[3]):
        return 0
    else: # x<y
        return -1

if len(sys.argv) == 3:
    infile = open(sys.argv[1])
    tot = sys.argv[2]
else:
    print "usage\nspaceOutSites.py <infile> <tot desired rows>"
    exit()

myset = [ line.split() for line in infile ]
    
myset.sort(cmp=my_compare1st)
    
setnum = 0
for s in range(1, int(tot)+1):
    if setnum < len(myset) and int(myset[setnum][0]) == s:
        print "\t".join(myset[setnum])
        setnum = setnum + 1
    else:
        print s
