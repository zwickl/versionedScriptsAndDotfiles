#!/usr/bin/env python
import os
import sys
import string
import StringIO

onlyCount = None
if len(sys.argv) > 2:
	filename = sys.argv[1]
	file = open(filename, 'rU')
else:
	print "usage\nfindAAColumns.py <filename> <AA to look for>"
        exit()

target = sys.argv[2]

for l in file:
        if 'XXX' in l:
	#if 1:
		num = 0
		for i in l[4:]:
			#print i
			num = num + 1
			if i == target:
				#print num
				print "incl", num*3-2, "-", num*3, ";"
				

#        if 'matrix' in l.lower():
#            ready = 1
#            print l
#            continue
#        if ready == 0:
#            print l
#        if ready and l.isspace() == 0:
#            #print l
#            if ';' in l:
#                print l
#                break
#            taxonNum = taxonNum + 1
#            baseCount = 0
#            lnum = lnum + 1
#            #print lnum
#            newline = ""
#            s=l.split()
#            #print s
#            name = s[0]
#            line = s[1]
#            pos = 0
#            while pos < len(line):
#                codon = line[pos:pos+3]
#                num = 0
#                if codon[0] in bases:
#                    num = num + 1
#                if codon[1] in bases:
#                    num = num + 1
#                if codon[2] in bases:
#                    num = num + 1
#                    
#                if num != 0 and num != 3:
#                    #print num, pos, codon
#                    newline += "???"
#                else:
#                    newline += codon
#                    baseCount = baseCount + num
#                    #print codon
#                #print pos
#                pos = pos+3
#            if baseCount > maxBaseCount:
#                maxBaseCount = baseCount
#            if onlyCount:
#                print name, "\t[", taxonNum, " ", baseCount, "]"
#            else:
#                print name, "\t[", taxonNum, " ", baseCount, "]\t", newline
#
#print "[Max nonambig bases = ", maxBaseCount, "]"
#
#if onlyCount == None:           
#    for l in file:
#        print l
#    
##print "done"
