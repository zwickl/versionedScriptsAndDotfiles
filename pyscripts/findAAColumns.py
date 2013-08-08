#!/usr/bin/env python
import sys

onlyCount = None
if len(sys.argv) > 2:
    filename = sys.argv[1]
else:
    print "usage\nfindAAColumns.py <filename> <AA to look for>"
    sys.exit()

target = sys.argv[2]

for l in open(filename, 'rU'):
    if 'XXX' in l:
        num = 0
        for i in l[4:]:
            num += 1
            if i == target:
                print "incl", num*3-2, "-", num*3, ";"
				
