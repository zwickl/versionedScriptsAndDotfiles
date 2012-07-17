#!/usr/bin/env python
import sys

onlyCount = None
if len(sys.argv) < 5:
    print "usage\nmakeCon.py <pos/neg (p or n)> <total taxa> <at least numbers to constrain ...>"
    exit(1)

tot = int(sys.argv[2])

conNums = []
for a in range(3, int(len(sys.argv))):
    conNums.append(int(sys.argv[a]))

pos = 0
if sys.argv[1] != 'p' and sys.argv[1] != 'n':
    print "second argument must be p or n"
else:
    if sys.argv[1] == 'p':
        print "+((",
    else:
        print "-((",
    
nums = range(1, tot+1)
toCon = []
for c in conNums:
    if c > tot:
        raise RuntimeError("taxon ", c, " is greater than the total number of taxa ", tot)
    toCon.append(c)
    nums.remove(c)

print toCon[0],
toCon.remove(toCon[0])

for i in toCon:
    print ", %d" % i,
print ")",
for i in nums:
    print ", %d" % i,
print ");"
