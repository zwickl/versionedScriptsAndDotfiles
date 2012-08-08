#!/usr/bin/env python
import sys
import re
import copy
import os
from argparse import ArgumentParser

from Bio import AlignIO
from Bio.Alphabet import IUPAC
from Bio.Align import MultipleSeqAlignment
from Bio.Nexus import Nexus

#use argparse module to parse commandline input
parser = ArgumentParser(description='extract sequences from a fasta file')

#add possible arguments
#flag
parser.add_argument('-i', '--interleave', dest='interleave', action='store_true', default=False,
                    help='interleave the nexus output matrix')

parser.add_argument('-q', '--quiet', dest='quiet', action='store_true', default=False,
                    help='output less crap to stderr')

parser.add_argument('-n', '--num-taxa', dest='numTax', type=int, default=None,
                    help='only consider alignments with the specified number of taxa')

#variable number of arguments
parser.add_argument('filenames', nargs='*', default=[], 
                    help='a list of filenames to search')



#now process the command line
options = parser.parse_args()

if len(options.filenames) == 0:
    raise RuntimeError("Enter alignment filenames to concatenate")
else:
    files = options.filenames

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
    
    if options.numTax is None or len(thisAlign) == options.numTax:
        endbase = startbase + thisAlign.get_alignment_length() - 1
        charsetString += "charset %s.%d = %d - %d;\n" % (re.sub('.*/', '', files[num-1]), num, startbase, endbase)
        startbase = endbase + 1
        num = num + 1
        
        allAlignments.append(thisAlign)

charsetString += "end;\n"

#the very simply old way, assuming equal # identically named seqs in all alignments
if oldConcat:
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
        if not options.quiet:
            sys.stderr.write("%d sequences in alignment\n" % len(thisDict))
        alignDicts.append((thisDict, align))
    sys.stderr.write("%d alignments read\n" % len(alignDicts))

    allNames = set()
    finalAlignSeqs = []
    for mydict in alignDicts:
        for seq in mydict[0].items():
            if seq[0] not in allNames:
                allNames |= set([seq[0]])
                finalAlignSeqs.append(copy.deepcopy(mydict[0][seq[0]]))
                finalAlignSeqs[-1]._set_seq('')
        #allNames |= set(mydict[0].keys())

    sys.stderr.write("%d names across all alignments\n" % len(allNames))

    for mydict in alignDicts:
        dummy = 'N' * mydict[1].get_alignment_length()
        for seq in finalAlignSeqs:
            if seq.name in mydict[0]:
                seq._set_seq(seq.seq + mydict[0][seq.name].seq)
            else:
                #seq.__add__ seems like it should work here, but I can't get it to
                seq._set_seq(seq.seq + dummy)

    finalAlign = MultipleSeqAlignment(finalAlignSeqs)

    if options.interleave:
        temp = '.temp.nex'
        AlignIO.write(finalAlign, temp, "nexus")
        backIn = Nexus.Nexus(temp)
        backIn.write_nexus_data(filename=sys.stdout, interleave=True)
        os.remove(temp)
    else:
        AlignIO.write(finalAlign, sys.stdout, "nexus")


sys.stdout.write("%s\n" % charsetString)

