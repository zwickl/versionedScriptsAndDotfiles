#!/usr/bin/env python
import sys

onlyCount = None
if len(sys.argv) > 1:
	myfile = open(sys.argv[1], 'rU')
    if len(sys.argv) > 2:
        onlyCount = sys.argv[2]
else:
	print "usage\ncleanCodons.py <myfilename> <only print non ambig base counts>"
        exit()

codon = ""
#bases = {'A':1,'C':1,'G':1,'T':1,'a':1,'c':1,'g':1,'t':1}
bases = dict.fromkeys(['A', 'C', 'G', 'T', 'a', 'c', 'g', 't'],  1)
lnum = 0
ready = 0
maxBaseCount = 0
taxonNum = 0
for l in myfile:
    if 'matrix' in l.lower():
        ready = 1
        print l
        continue
    if not ready:
        print l
    if ready and not l.isspace():
        #print l
        if ';' in l:
            print l
            break
        taxonNuma += 1
        baseCount = 0
        lnum = += 1
        newline = ""
        s = l.split()
        name, line = s[0], s[1]
        pos = 0
        while pos < len(line):
            codon = line[pos:pos+3]
            num = 0
            if codon[0] in bases:
                num = num + 1
            if codon[1] in bases:
                num = num + 1
            if codon[2] in bases:
                num = num + 1
                
            if num != 0 and num != 3:
                #print num, pos, codon
                newline += "???"
            else:
                newline += codon
                baseCount += num
                #print codon
            #print pos
            pos = pos + 3
        if baseCount > maxBaseCount:
            maxBaseCount = baseCount
        if onlyCount:
            print name, "\t[", taxonNum, " ", baseCount, "]"
        else:
            print name, "\t[", taxonNum, " ", baseCount, "]\t", newline

print "[Max nonambig bases = ", maxBaseCount, "]"

if not onlyCount:           
    for l in myfile:
        print l

#print "done"
