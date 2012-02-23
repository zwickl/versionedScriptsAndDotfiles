#!/usr/bin/env python

import sys
import dendropy

sys.stderr.write('Reading %s ...\n' % sys.argv[1])

intree= dendropy.Tree.get_from_path(sys.argv[1], "nexus")

if len(sys.argv) > 2:
    out = open(sys.argv[2], 'w')
else:
    out = sys.stdout

intree.write(out, "newick")

