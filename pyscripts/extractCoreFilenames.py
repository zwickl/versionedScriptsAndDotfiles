#!/usr/bin/env python

import sys
from dzutils import extract_core_filename

if len(sys.argv) > 1:
    inf = open(sys.argv[1], 'rb')
else:
    inf = sys.stdin

for line in inf:
    #print extract_core_filename(line)
    print extract_core_filename(line, no_nchar=True)


