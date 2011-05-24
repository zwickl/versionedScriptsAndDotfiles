#!/usr/bin/env python
import os
import sys
import string
import StringIO


if len(sys.argv) < 3:
    print "usage\nsumNucSitelikes.py <# sites per codon (2 or 3)> <filename>"
    exit()
else:
    if sys.argv[1] == '2' or sys.argv[1] == '3':
        sitesPerCodon = int(sys.argv[1])   
    else:
        print "usage\nsumNucSitelikes.py <# sites per codon (2 or 3)> <filename>"
        exit()   
    filename = sys.argv[2]
    file = open(filename, 'rU')

poo = ""
lnLs = []
print "codon#   lnL"
for i in file:
    poo = i
    spl = poo.split()
    if spl[0].isdigit() > 0:
        if spl[1] == "gap":
            lnLs.append(0.0)
        else:
            lnLs.append(float(spl[1]))
            
        if len(lnLs) == sitesPerCodon:
            if sitesPerCodon == 3:
                sum = lnLs[0] + lnLs[1] + lnLs[2]
                codNum = int(spl[0]) / 3
                print  codNum, "\t", sum
            else:
                sum = lnLs[0] + lnLs[1]
                codNum = int(spl[0]) / 2
                print  codNum, "\t", sum
            lnLs = []
    