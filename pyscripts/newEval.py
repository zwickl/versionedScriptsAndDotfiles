#!/usr/bin/env python
import os
import string

def islogfile(x): 
    return x.endswith('.log00.log')

def numeric_compare(x, y):
    f = x[:2]
    s = y[:2]
    return int(f) - int(s)

runnum = 1
files = os.listdir('.')
#files = filter(islogfile, files)
files = [ f for f in files if islogfile(f) ]
files.sort(cmp=numeric_compare)
if not files:
    print "Sorry, no log files found"

expected = [ line for line in open('expected.scr', 'rb') ]

num = 0
for i in files:
    if num >= len(expected):
        break
    if i[0].isdigit():
        scores = []
        infile = open(i, 'r')
        line = infile.readline()
        while line != "":
            if line.startswith('Final'):
                divided = line.split()
                scores.append(string.atof(divided[1]))
            line = infile.readline()
#        print scores
    #    sum = string.atof(expected[num]) + string.atof(divided[1])
        if scores:
            mySum = string.atof(expected[num]) + max(scores)
            #print '%(f)s %(e)f %(g)s' % {'f': i, "e": string.atof(expected[num]), "g": divided[1]}
            if abs(mySum) < 0.05:
                print '%(f)25s %(s)15f PASS ' % {'f': i, "s": mySum}, 
            else:
                print '%(f)25s %(s)15f FAIL' % {'f': i, "s": mySum},
            print expected[num], scores
        else:
            print i, "no score found in infile?"
        infile.close()
        num = num + 1






