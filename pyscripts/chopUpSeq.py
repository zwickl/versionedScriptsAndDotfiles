#!/usr/bin/env python
import os
import sys
import string
import StringIO

offset = 0
if len(sys.argv) > 2:
    filename = sys.argv[1]
    file = open(filename, 'rU')
    perLine = sys.argv[2]
    if len(sys.argv) > 3:
        offset = int(sys.argv[3])
else:
    print "usage\nchopUpSeq.y <filename> <#per line> <optional codon NUMBER offset for first line>"
    exit()

seq = file.read()
    
joined = seq.replace("\n", "")
#print joined

i = 0
pl = int(perLine)

if offset:
    print joined[0:offset-1]

while(i * pl + offset - 1< len(joined)):
    if i * pl + pl + offset > len(joined):
        print joined[i * pl - 1 + offset:len(joined)]
    else:
        print joined[i * pl - 1 + offset:i * pl + pl - 1 + offset]
    i = i + 1
    
#print seq

