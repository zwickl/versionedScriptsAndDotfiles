#!/usr/bin/env python
import sys

if len(sys.argv) != 2:
    print "usage\nsitelikeTable.py <sitelike infilename>"
    exit()
else:
    infile = open(sys.argv[1], 'rb')

trees = []
totals = []
thistree = []
for line in infile.readlines():
    if not 'Tree' in line:
        #depending on the fact that sitelike lines start with tabs and total tree like lines don't
        if line[0] == '\t':
            thistree.append(line.split())
        else:
            totals.append(line.split())
            trees.append(thistree)
            thistree = []

if len(thistree):
    trees.append(thistree)

numStr = '\t'.join([ "%s" % n[0] for n in totals ])
scoreStr = '\t'.join([ "%f" % float(s[1]) for s in totals ])

print "tree#", numStr
print "scores", scoreStr

for site in xrange(len(trees[0])):
    print site, '\t'.join([tree[site][1] for tree in trees])

