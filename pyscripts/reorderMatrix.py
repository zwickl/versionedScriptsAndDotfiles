#!/usr/bin/env python
import sys
import random

rows = []

elements = [["A", "Asp"]]

#elements.append(["A","Ala"])
elements.append(["C","Cys"])
elements.append(["D","Asp"])
elements.append(["E","Glu"])
elements.append(["F","Phe"])
elements.append(["G","Gly"])
elements.append(["H","His"])
elements.append(["I","Ile"])
elements.append(["K","Lys"])
elements.append(["L","Leu"])
elements.append(["M","Met"])
elements.append(["N","Asn"])
elements.append(["P","Pro"])
elements.append(["Q","Gln"])
elements.append(["R","Arg"])
elements.append(["S","Ser"])
elements.append(["T","Thr"])
elements.append(["V","Val"])
elements.append(["W","Trp"])
elements.append(["Y","Tyr"])




rows.append("A", "Asp", [)















n = float(sys.argv[1])
chars = range(int(n/3))
num_sets = int(sys.argv[2])
prop = float(sys.argv[3])
a = int(n*prop/3.0)
#print prop, a
for i in range(num_sets):
    s = random.sample(chars, a)
    B =  " ".join(["%d %d %d" %(3*i+1, 3*i+2, 3*i+3) for i in s])
    print "Exset * bogus = %s;" % B
    
