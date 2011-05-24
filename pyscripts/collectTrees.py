#!/usr/bin/env python
import os
import sys
import glob
import string
import StringIO

if len(sys.argv) > 1:
    if len(sys.argv) > 2:
        #assume that an actual list of files was passed in, as would be if a glob were expanded
        #by the shell, i.e., collectTrees.py *.tre
        #print "using multiple cmd line arguments (or shell glob) as list of tree files ..."
        sys.stderr.write('using multiple cmd line arguments (or shell glob) as list of tree files ...\n')
        files = sys.argv[1:]
    else:
        #try interpretting as a literal glob pattern being passed in , i.e., collectTrees.py "*.tre"
        #print "using cmd line argument as literal glob ..."
        sys.stderr.write('using cmd line argument as literal glob ...\n')
        patt = sys.argv[1]
        files = glob.glob(patt)
        #print files
        if len(files) < 2:
            #print "nope.  using cmd line argument as path or partial filename pattern for .tre files"
            sys.stderr.write('nope.  using cmd line argument as path or partial filename pattern for .tre files\n')
            #if that didn't return anything besides the patt itself (as expected it nothing matches)
            #try interpretting as a path at which to look for .tre files
            patt = sys.argv[1] + "*.tre"
            files = glob.glob(patt)
            if len(files) < 2:
                #print "at least 2 tree files must exist"
                sys.stderr.write('at least 2 tree files must exist\n')
                exit()
else:
    #print "as command line arguments, enter\n\ta list of tree files\n\ta literal glob surrounded by quotes\n\tpath or partial pattern to which *.tre can be appended"
    sys.stderr.write('as command line arguments, enter\n\ta list of tree files\n\ta literal glob surrounded by quotes\n\tpath or partial pattern to which *.tre can be appended\n')
    exit()

file1 = open(files[0], "rU")

#get the translate block from the first file
for line in file1:
    if line.find("tree ") < 0:
        print line,
    else:
        break

for f in files:
    eachFile = open(f, "rU")
    for indLine in eachFile:
        if indLine.lower().find("tree") > -1 and indLine.lower().find("=") > -1:
            print indLine,

print "end;"

