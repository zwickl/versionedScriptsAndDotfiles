#!/usr/bin/env python
from dzutils import *
import sys
import re
import itertools
from math import *
import argparse

#use argparse module to parse commandline input
parser = argparse.ArgumentParser(description='extract sequences from a fasta file')

#add possible arguments
parser.add_argument('-v', '--invert-match', dest='invertMatch', action='store_true', default=False,
                    help='invert the sense of the match (default false)')

parser.add_argument('--single-copy', dest='onlySingleCopy', action='store_true', default=False,
                    help='only include single copy clusters (default false)')

parser.add_argument('-s', '--sort', dest='sortOutput', action='store_true', default=False,
                    help='alphanumerically sort the sequences by name (default false)')

parser.add_argument('-f', '--patternfile', dest='patternFile', type=str, default=None, 
                    help='file from which to read patterns (you must still pass a pattern on the command line, which is ignored)')

parser.add_argument('pattern',
                    help='a quoted regular expression to search sequence names for')

parser.add_argument('--range', dest='baseRange', nargs=2, type=int, default=[1, 9999999], metavar=('startbase', 'endbase'),
                    help='range of cluster sizes (number of members)')

parser.add_argument('filenames', nargs='*', default=[], 
                    help='a list of filenames to search (none for stdin)')

#now process the command line
parsed = parser.parse_args()

clust = parse_blink_output(parsed.filenames[0])

if parsed.patternFile is not None:
    sys.stderr.write('reading patterns from file %s ...\n' % parsed.patternFile)
    pf = open(parsed.patternFile, 'rb')
    taxPatterns = [ line.strip() for line in pf ]
    sys.stderr.write('patterns: %s\n' % str(taxPatterns))
else:
    taxPatterns = [parsed.pattern]
seqFiles = parsed.filenames

if parsed.invertMatch:
    matchedRecs = set(clust)
else:
    matchedRecs = set()

for cpat in taxPatterns:
    for c in clust:
        if len(c) >= parsed.baseRange[0] and len(c) <= parsed.baseRange[1]:
            if parsed.onlySingleCopy is False or c.is_single_copy():
                match = c.contains_matching_taxon(cpat)
                if match:
                    if parsed.invertMatch:
                        if c in matchedRecs:
                            matchedRecs.remove(c)
                    else:
                        matchedRecs.add(c)
         
for r in matchedRecs:
     r.output()

'''

for c in clust:
    found = c.contains_matching_taxon(parsed.pattern)
    if len(c) >= parsed.baseRange[0] and len(c) <= parsed.baseRange[1]:
        if found and parsed.invertMatch is False:
            c.output()
        elif found == False and parsed.invertMatch:
            c.output()

'''
