#!/usr/bin/env python
import sys

if len(sys.argv) < 3:
    print "usage\nsumNucSitelikes.py <# sites per codon (2 or 3)> <filename>"
    sys.exit()
else:
    if sys.argv[1] == '2' or sys.argv[1] == '3':
        sitesPerCodon = int(sys.argv[1])   
    else:
        print "usage\nsumNucSitelikes.py <# sites per codon (2 or 3)> <filename>"
        sys.exit()   
    filename = sys.argv[2]

poo = ""
lnLs = []
print "codon#   lnL"
for i in open(filename, 'rU'):
    poo = i
    spl = poo.split()
    if spl[0].isdigit() > 0:
        if spl[1] == "gap":
            lnLs.append(0.0)
        else:
            lnLs.append(float(spl[1]))
            
        if len(lnLs) == sitesPerCodon:
            if sitesPerCodon == 3:
                mySum = lnLs[0] + lnLs[1] + lnLs[2]
                codNum = int(spl[0]) / 3
                print  codNum, "\t", mySum
            else:
                mySum = lnLs[0] + lnLs[1]
                codNum = int(spl[0]) / 2
                print  codNum, "\t", mySum
            lnLs = []
    
