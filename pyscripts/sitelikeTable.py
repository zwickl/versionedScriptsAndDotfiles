#!/usr/bin/env python
import sys

if len(sys.argv) != 2:
    print "usage\nsitelikeTable.py <sitelike filename>"
    exit()
else:
    filename = sys.argv[1]
    file = open(filename, 'rU')

#get rid of the first header line    
file.readline()
trees = []
totals = []
thistree = []
for line in file.readlines():
    #depending on the fact that sitelike lines start with tabs and total tree like lines don't
    if line[0] == '\t':
        thistree.append(line.split())
    else:
        totals.append(line.split())
        trees.append(thistree)
        thistree = []
#        print line

if len(thistree) > 0:
    trees.append(thistree)

print "tree#\t",

for t in range(0, len(trees)):
    print "%d\t" % int(totals[t][0]) ,
print "\nscore\t",
for t in range(0, len(trees)):
    print "%f\t" % float(totals[t][1]) , 
print "\ntree\t",
for t in range(0, len(trees)):
    print "%s\t" % totals[t][2] , 
print

for site in range(0, len(trees[0])):
    print "%d\t" % int(trees[t][site][0]),
    for t in range(0, len(trees)):
        print "%f\t" % float(trees[t][site][1]),
    print
