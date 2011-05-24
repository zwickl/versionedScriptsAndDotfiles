#!/usr/bin/env python
import os
import sys
import subprocess
import string
from StringIO import StringIO
import re
import itertools
from optparse import OptionParser

from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIStandalone

parser = OptionParser(add_help_option=False)

parser.add_option("-h", "--help", dest="helpflag")
parser.add_option("-v", "--invert-match", action="store_true", default=False, dest="invertFlag", help="invert match sense")
#parser.add_option("-i", "--inputfile", dest="filename", type="string", help="input file", metavar="FILE")
#parser.add_option("-r", "--rows", dest="nrows", type="int", help="number of rows per chunk", metavar="#")
#parser.add_option("-f", "--firstcols", dest="firstCols", type="string", action="callback", callback=comma_split, help="columns to use in first row set", metavar="# # ...")
#parser.add_option("-c", "--cols", dest="cols", type="string", action="callback", callback=comma_split, help="columns to use in successive row sets", metavar="# # ...")
#parser.add_option("-t", "--transpose", action="store_true", dest="transflag", default=False, help="transpose matrix")

(options, args) = parser.parse_args()

'''
if options.helpflag != None:
    print_usage("")

if options.filename == None:
print_usage("you must pass the -i option")
if options.nrows == None:
print_usage("you must pass the -r option")
if options.firstCols == None:
print_usage("you must pass the -f option")
if options.cols == None:
print_usage("you must pass the -c option")
'''

invertMatch = options.invertFlag

if len(sys.argv) < 3:
    sys.stderr.write("Enter a quoted sequence name regex and sequence file(s)\n")
    exit(0)
else:
    seqFiles = args[1:]
    seqPattern = args[0]

sys.stderr.write("matching pattern %s ...\n" % seqPattern)
try:
    pat = re.compile(seqPattern)
except:
    sys.stderr.write("problem compiling regex pattern %s\n" % seqPattern)
    exit(1)

sys.stderr.write("Parsing sequence files %s ...\n" % str(seqFiles))

for oneSeqFile in seqFiles:
    try:
        allRecs = [ rec for rec in SeqIO.parse(oneSeqFile, "fasta", alphabet=IUPAC.ambiguous_dna) ]
    except IOError:
        sys.stderr("error reading file %s!" % oneSeqFile)
        exit(1)

    matchedSeqs = []
    for seq in allRecs:
        match = pat.search(seq.description)
        if ( match is not None and invertMatch is False ) or ( match is None and invertMatch is True ):
        #if match is not None:
            matchedSeqs.append(seq)
    
    if len(matchedSeqs) > 0:
        sys.stderr.write("matched %d sequences in %s\n" % (len(matchedSeqs), oneSeqFile))
    else:
        sys.stderr.write("no sequences matched in %s!\n" % oneSeqFile)

    SeqIO.write(matchedSeqs, sys.stdout, "fasta")

'''
if len(sys.argv) < 3:
    sys.stderr.write("Enter a sequence filename and a quoted sequence name regex(s)\n")
    exit(0)
else:
    allSeqFile = sys.argv[1]
    seqPatterns = sys.argv[2:]

sys.stderr.write("Parsing sequence file %s ...\n" % allSeqFile)
allRecs = [ rec for rec in SeqIO.parse(allSeqFile, "fasta", alphabet=IUPAC.ambiguous_dna) ]

matchedSeqs = []
for patStr in seqPatterns:
    sys.stderr.write("matching pattern %s ...\n" % patStr)
    try:
        pat = re.compile(patStr)
    except:
        sys.stderr.write("problem compiling regex pattern %s" % patStr)
        exit(1)

    patMatches = []
    for seq in allRecs:
        match = pat.search(seq.description)
        if match is not None:
            patMatches.append(seq)
            sys.stderr.write("matched %d patterns\n" % len(patMatches))
            matchedSeqs.extend(patMatches)
            if len(matchedSeqs) > 0:
                sys.stderr.write("matched sequences: %s\n" % ", ".join([ str(seq.name) for seq in matchedSeqs ]))
            else:
                sys.stderr.write("no sequences matched!\n")
                exit(0)

sys.stderr.write("outputting %s sequences\n" % len(matchedSeqs))
SeqIO.write(matchedSeqs, sys.stdout, "fasta")
'''
