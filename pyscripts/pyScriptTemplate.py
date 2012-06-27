#!/usr/bin/env python
import sys
from argparse import ArgumentParser

#use argparse module to parse commandline input
parser = ArgumentParser(description='extract sequences from a fasta file')

#use mutExclGroup.add_argument instead of parser.add_argument for mutually exclusive options
#mutExclGroup = parser.add_mutually_exclusive_group()

#add possible arguments
#flag
parser.add_argument('-v', '--invert-match', dest='invertMatch', action='store_true', default=False,
                    help='invert the sense of the match (default false)')

#string
parser.add_argument('-f', '--patternfile', dest='patternFile', type=str, default=None, 
                    help='file from which to read patterns (you must still pass a pattern on the command line, which is ignored)')

#multiple arguments
parser.add_argument('--range', dest='baseRange', nargs=2, type=int, default=[1, 9999999], metavar=('startbase', 'endbase'),
                    help='range of cluster sizes (number of members)')

#single number value
parser.add_argument('-mp', '--min-match-prop', dest='minMatchProportion', type=float, default=0.0,
                    help='proportion of hit that must overlap query (default 0.0)')

#variable number of arguments
parser.add_argument('filenames', nargs='*', default=[], 
                    help='a list of filenames to search (none for stdin)')

#now process the command line
options = parser.parse_args()

