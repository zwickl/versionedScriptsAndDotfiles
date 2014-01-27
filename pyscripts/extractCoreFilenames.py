#!/usr/bin/env python

import sys
from argparse import ArgumentParser
from dzutils import extract_core_filename

parser = ArgumentParser()

parser.add_argument('filename', nargs='*', type=str, default=None)

parser.add_argument('-c', '--with-chars', action='store_true', default=False, help='include the part of the core name including the number of characters')

args = parser.parse_args()

if args.filename:
    inf = open(args.filename, 'rb')
else:
    inf = sys.stdin

for line in inf:
    sys.stdout.write('%s\n' % extract_core_filename(line, no_nchar=not args.with_chars))


