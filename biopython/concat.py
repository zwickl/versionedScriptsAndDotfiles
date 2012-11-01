#!/usr/bin/env python
import sys
import re
import copy
import os
from argparse import ArgumentParser

from Bio import AlignIO
from Bio import Seq
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

if not options.filenames:
    raise RuntimeError("Enter alignment filenames to concatenate")

oldConcat = False

allAlignments = []

charsetString = "begin sets;\n"
charpartString = "charpartition concat = "
num, startbase, endbase = 1, 1, 1

for filename in options.filenames:
    with open(filename, "rb") as afile:
        thisAlign = AlignIO.read(afile, "nexus", alphabet=IUPAC.ambiguous_dna)
    
    if not options.numTax or len(thisAlign) == options.numTax:
        endbase = startbase + thisAlign.get_alignment_length() - 1
        #charsetString += "charset %s.%d = %d - %d;\n" % (re.sub('.*/', '', files[num-1]), num, startbase, endbase)
        charsetString += "charset c%d = %d - %d; [%s]\n" % (num, startbase, endbase, re.sub('.*/', '', options.filenames[num-1]))
        charpartString += "%d:c%d, " % (num, num)
        startbase = endbase + 1
        num = num + 1
        
        allAlignments.append(thisAlign)

charsetString += charpartString 
charsetString += ";\nend;\n"

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

    #figure out all necessary taxa in the final alignment
    allNames = set()
    for mydict, alignment in alignDicts:
        allNames |= set(mydict.iterkeys())

    '''
    finalAlignSeqs = []
    for mydict, alignment in alignDicts:
        for name, seq in mydict.items():
            if name not in allNames:
                allNames.add(name)
                finalAlignSeqs.append(copy.deepcopy(seq))
                finalAlignSeqs[-1]._set_seq('')
    '''

    sys.stderr.write("%d names across all alignments\n" % len(allNames))

    #collect the sequences for each taxa for each alignment, or a dummy string of N's if they are missing from an alignment
    #rawStrings = dict.fromkeys(allNames)
    rawStrings = {}
    
    #loop over alignments
    for mydict, alignment in alignDicts:
        dummy = 'N' * alignment.get_alignment_length()
        #loop over the list of sequences that we're collecting
        for name in allNames:
            if name in mydict:
                rawStrings.setdefault(name, []).append(str(mydict[name].seq))
            else:
                rawStrings.setdefault(name, []).append(dummy)

    #this is a little cheesy, but find the first alignment that contains each taxon and use the 
    #corresponding SeqRecord as a template to paste the new full alignments into
    finalAlignSeqs = []
    for name in sorted(allNames):
        for mydict, alignment in alignDicts:
            if name in mydict:
                finalSeq = copy.deepcopy(mydict[name])
                finalSeq._set_seq(Seq.Seq(''.join(rawStrings[name]), finalSeq.seq.alphabet))
                finalAlignSeqs.append(finalSeq)
                break

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

