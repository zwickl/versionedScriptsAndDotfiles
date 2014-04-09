#!/usr/bin/env python
from dzutils import parse_blink_output
import sys
#from math import *
import argparse

#use argparse module to parse commandline input
parser = argparse.ArgumentParser(description='extract sequences from a fasta file')

#add possible arguments
parser.add_argument('-v', '--invert-match', action='store_true', default=False,
                    help='invert the sense of the match (default false)')

parser.add_argument('--single-copy', action='store_true', default=False,
                    help='only include single copy clusters (default false)')

parser.add_argument('-s', '--sort', action='store_true', default=False,
                    help='sort clusters by cluster number (default false)')

parser.add_argument('--list-only', action='store_true', default=False,
                    help='only output the numbers of the clusters that match, not the actual sequences in the clusters')

parser.add_argument('-f', '--pattern-file', type=str, default=None, 
                    help='file from which to read patterns (you must still pass a pattern on the command line, which is ignored)')

parser.add_argument('pattern', type=str, 
                    help='a quoted regular expression to search sequence names or cluster numbers for')

parser.add_argument('-n', '--nums', dest='numbers', action='store_true', default=False,
                    help='assume that the patterns refer to the cluster numbers, not taxon names (default false)')

parser.add_argument('--range', nargs=2, type=int, default=[1, 9999999], metavar=('smallest', 'largest'),
                    help='range of cluster sizes (number of members)')

parser.add_argument('filenames', nargs='*', default=[], 
                    help='a list of filenames to search (none for stdin)')

#now process the command line
options = parser.parse_args()

clust = parse_blink_output(options.filenames[0])

if options.pattern_file:
    sys.stderr.write('reading patterns from file %s ...\n' % options.pattern_file)
    taxPatterns = [ line.strip() for line in open(options.patternFile, 'rb') ]
    sys.stderr.write('patterns: %s\n' % str(taxPatterns))
else:
    taxPatterns = [options.pattern]

matchedRecs = set(clust) if options.invert_match else set()

for cpat in taxPatterns:
    for c in clust:
        if options.range[0] <= len(c) <= options.range[1]:
            if not options.single_copy or c.is_single_copy():
                if not options.numbers:
                    match = c.contains_matching_taxon(cpat)
                    if match:
                        if options.invert_match:
                            if c in matchedRecs:
                                matchedRecs.remove(c)
                        else:
                            matchedRecs.add(c)
                else:
                    if c.number == int(cpat):
                        if options.invert_match:
                            if c in matchedRecs:
                                matchedRecs.remove(c)
                        else:
                            matchedRecs.add(c)

if options.sort:
    matchedRecs = sorted(list(matchedRecs), key=lambda n:n.number)

for r in matchedRecs:
    if options.list_only:
        sys.stdout.write('%s\n' % r.number)
    else:
        r.output()

'''

for c in clust:
    found = c.contains_matching_taxon(options.pattern)
    if len(c) >= options.baseRange[0] and len(c) <= options.baseRange[1]:
        if found and not options.invert_match:
            c.output()
        elif not found and options.invert_match:
            c.output()

'''
