#!/usr/bin/env python
import os
import sys
import subprocess
import string
from StringIO import StringIO
import re
import itertools
import copy
import time

from Bio import AlignIO
from Bio.Alphabet import IUPAC

if len(sys.argv) < 2:
    sys.stderr.write("Enter alignment filenames to concatenate\n")
    exit(0)
else:
    files = sys.argv[1:]

allAlignments = []

for filename in files:
    try:
        afile = open(filename, "rb")
    except IOError:
        exit("could not read file %s" % filename)
    allAlignments.append(AlignIO.read(afile, "nexus", alphabet=IUPAC.ambiguous_dna))

num = 1
concat = allAlignments[0]
startbase = 1
endbase = allAlignments[0].get_alignment_length()
charsetString = "begin assumptions;\n"
charsetString += "charset gene%d = %d - %d;\n" % (num, startbase, endbase)
num = num + 1
startbase = endbase + 1
for align in allAlignments[1:]:
    concat += align
    endbase = startbase + align.get_alignment_length() - 1
    charsetString += "charset gene%d = %d - %d;\n" % (num, startbase, endbase)
    startbase = endbase + 1
    num = num + 1
charsetString += "end;\n"

AlignIO.write(concat, sys.stdout, "nexus")

print charsetString
