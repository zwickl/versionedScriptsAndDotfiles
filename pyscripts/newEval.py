#!/usr/bin/env python
import os
import string

def numeric_compare(x, y):
    f = x[:2]
    s = y[:2]
    return int(f) - int(s)

runnum = 1
files = os.listdir('.')
files = [ f for f in files if f.endswith('.log00.log') ]
files.sort(cmp=numeric_compare)
if not files:
    print "Sorry, no log files found"

expected = [ line for line in open('expected.scr', 'rb') ]

num = 0
for i in files:
    while num < len(expected):
    if i[0].isdigit():
        scores = []
        with open(i, 'r') as infile:
            line = infile.readline()
            while line != "":
                if line.startswith('Final'):
                    divided = line.split()
                    scores.append(string.atof(divided[1]))
                line = infile.readline()
            if scores:
                mySum = string.atof(expected[num]) + max(scores)
                if abs(mySum) < 0.05:
                    print '%(f)25s %(s)15f PASS ' % {'f': i, "s": mySum}, 
                else:
                    print '%(f)25s %(s)15f FAIL' % {'f': i, "s": mySum},
                print expected[num], scores
            else:
                print i, "no score found in infile?"
            num = num + 1






