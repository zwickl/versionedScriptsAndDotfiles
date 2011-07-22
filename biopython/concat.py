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
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Align import MultipleSeqAlignment

if len(sys.argv) < 2:
    sys.stderr.write("Enter alignment filenames to concatenate\n")
    exit(0)
else:
    files = sys.argv[1:]

oldConcat = False

allAlignments = []

charsetString = "begin assumptions;\n"
num = 1
startbase = 1
endbase = 1

for filename in files:
    try:
        afile = open(filename, "rb")
    except IOError:
        exit("could not read file %s" % filename)
    thisAlign = AlignIO.read(afile, "nexus", alphabet=IUPAC.ambiguous_dna)
    
    endbase = startbase + thisAlign.get_alignment_length() - 1
    charsetString += "charset %s.%d = %d - %d;\n" % (re.sub('.*/', '', sys.argv[num]), num, startbase, endbase)
    startbase = endbase + 1
    num = num + 1
    
    allAlignments.append(thisAlign)

charsetString += "end;\n"

#the very simply old way, assuming equal # identically named seqs in all alignments
if oldConcat is True:
    concat = allAlignments[0]
    for align in allAlignments[1:]:
        concat += align

    AlignIO.write(concat, sys.stdout, "nexus")

else:
    #actually, this is a list of (name:seq dict, alignment) tuples
    alignDicts = []
    for align in allAlignments:
        thisDict = {}
        for seq in align:
            name = seq.name
            if name in thisDict:
                sys.stderr.write("sequence name %s already found in alignment:\n" % name)
                sys.stderr.write(align)
                exit(1)
            try:
                thisDict[name] = seq
            except KeyError:
                sys.stderr.write("problem adding sequence %s to alignment dictionary\n" % name)
                sys.stderr.write(align)
                exit(1)
        sys.stderr.write("%d sequences in alignment\n" % len(thisDict))
        alignDicts.append((thisDict, align))
    sys.stderr.write("%d alignments read\n" % len(alignDicts))

    allNames = set()
    finalAlignSeqs = []
    for dict in alignDicts:
        for seq in dict[0].items():
            if seq[0] not in allNames:
                allNames |= set([seq[0]])
                finalAlignSeqs.append(copy.deepcopy(dict[0][seq[0]]))
                finalAlignSeqs[-1]._set_seq('')
        #allNames |= set(dict[0].keys())

    sys.stderr.write("%d names across all alignments\n" % len(allNames))

    for dict in alignDicts:
        dummy = ''
        for i in range(0,dict[1].get_alignment_length()):
            dummy += 'N'
        for seq in finalAlignSeqs:
            name = seq.name
            try:
                toAppend = dict[0][name].seq
            except KeyError:
                toAppend = dummy
            seq._set_seq(seq.seq + toAppend)
            
            #seq.__add__(thisSeq)

    finalAlign = MultipleSeqAlignment(finalAlignSeqs)
    AlignIO.write(finalAlign, sys.stdout, "nexus")

sys.stdout.write("%s\n" % charsetString)

