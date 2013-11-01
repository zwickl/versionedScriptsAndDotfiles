#!/usr/bin/env python
import sys

if len(sys.argv) < 5:
    sys.exit("usage\nmakeCon.py <pos/neg (p or n)> <total taxa> <at least numbers to constrain ...>")

tot = int(sys.argv[2])

conNums = [ int(arg) for arg in sys.argv[3:] ]

if sys.argv[1] not in ['p', 'n']:
    print "second argument must be p or n"
else:
    if sys.argv[1] == 'p':
        print "+((",
    else:
        print "-((",
    
allNums = set(range(1, tot+1))
toCon = set(conNums)
if max(toCon) > max(allNums):
    sys.exit('numbers to constraint > than total number')

allNums -= toCon

print ", ".join([str(n) for n in sorted(toCon)]) + "), " + ", ".join([str(n) for n in sorted(allNums)]) + ");"

