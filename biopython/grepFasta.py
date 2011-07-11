#!/usr/bin/env python
import os
import sys
import subprocess
import string
from StringIO import StringIO
import re
import itertools
import argparse

from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIStandalone

#use argparse module to parse commandline input
parser = argparse.ArgumentParser(description='extract sequences from a fasta file')

#add possible arguments
parser.add_argument('-v', '--invert-match', dest='invertMatch', action='store_true', default=False,
                    help='invert the sense of the match (default false)')

parser.add_argument('--range', dest='baseRange', nargs=2, type=int, default=[1, -1], metavar=('startbase', 'endbase'),
                    help='range of alignment positions to output, start at 1, last position included, -1 for end')

parser.add_argument('pattern',
                    help='a quoted regular expression to search sequence names for')

parser.add_argument('filenames', nargs='*', default=[], 
                    help='a list of filenames to search (none for stdin)')

#now process the command line
parsed = parser.parse_args()

invertMatch = parsed.invertMatch
startBase = parsed.baseRange[0] - 1
endBase = parsed.baseRange[1]
print startBase, endBase
seqPattern = parsed.pattern
seqFiles = parsed.filenames

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
            matchedSeqs.append(seq[startBase:endBase])
    
    if len(matchedSeqs) > 0:
        sys.stderr.write("matched %d sequences in %s\n" % (len(matchedSeqs), oneSeqFile))
    else:
        sys.stderr.write("no sequences matched in %s!\n" % oneSeqFile)

    SeqIO.write(matchedSeqs, sys.stdout, "fasta")

