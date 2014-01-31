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

def extract_name_with_regex(pattern, name):
    match = re.search(pattern, name)
    if not match:
        sys.exit('pattern doesn\'t match name %s' % name)
    #multiple groups can be specified in the pattern, and if some don't match they will be None
    return ''.join([g for g in match.groups() if g])


#use argparse module to parse commandline input
parser = ArgumentParser(description='concatenate a number of alignments, matching up taxon names across them')

#add possible arguments
#flag
parser.add_argument('-i', '--interleave', action='store_true', default=False,
                    help='interleave the nexus output matrix')

parser.add_argument('-q', '--quiet', action='store_true', default=False,
                    help='output less crap to stderr')

parser.add_argument('--no-char-partition', action='store_true', default=False,
                    help='don\'t include a sets blocks with a charset for each gene and a charpartition')

parser.add_argument('--ignore-multiple-copy', action='store_true', default=False,
                    help='just skip over an alignment if multiple names map to the same matchable name (default is to error out)')

parser.add_argument('-n', '--num-taxa', type=int, default=None,
                    help='only consider alignments with the specified number of taxa')

parser.add_argument('--name-pattern', type=str, default='(.*)', 
                    help='regex to use to pull names to match across alignments from names within alignments. \
                            Uses capture groups, and if multiple they will be concatenated to make final names')

parser.add_argument('-r', '--taxa-range', nargs=2, default=[1, 9999999999], 
                    help='only consider alignments with between the specified range of taxa')

#variable number of arguments
parser.add_argument('filenames', nargs='*', default=[], 
                    help='a list of filenames to search')

#now use the tkinter gui, which may not work, or process the command line
use_tk_gui = False
if use_tk_gui:
    from tkinterutils import *
    root = Tk()
    gui = ArgparseGui(root, parser, height=1000, width=1600)
    
    root.wait_window(gui.frame)
    if gui.cancelled:
        sys.exit('cancelled ...')
    args = gui.make_commandline_list()
    print args

    options = parser.parse_args(args)
    print options
else:
    options = parser.parse_args()

if not options.filenames:
    raise RuntimeError("Enter alignment filenames to concatenate")

if options.num_taxa:
    minTax = maxTax = options.num_taxa
else:
    (minTax, maxTax) = options.taxa_range

alignDicts = []

charsetString = "begin sets;\n"
charpartString = "charpartition concat = "
num, startbase, endbase = 1, 1, 1

for filename in options.filenames:
    with open(filename, "rb") as afile:
        try:
            thisAlign = AlignIO.read(afile, "nexus", alphabet=IUPAC.ambiguous_dna)
        except ValueError:
            sys.stderr.write('problem reading alignment %s' % filename)
            raise
    
    thisDict = {}
    if minTax <= len(thisAlign) <= maxTax:
        for seq in thisAlign:
            name =  extract_name_with_regex(options.name_pattern, seq.name)
            if name in thisDict:
                sys.stderr.write("sequence name %s already found in alignment:\n" % name)
                sys.stderr.write("filename %s\n" % filename)
                if options.ignore_multiple_copy:
                    thisDict = {}
                    break
                else:
                    sys.exit(1)
            try:
                thisDict[name] = seq
            except KeyError:
                sys.stderr.write("problem adding sequence %s to alignment dictionary\n" % name)
                sys.stderr.write("filename %s\n" % filename)
                sys.exit(1)
        if thisDict:
            if not options.quiet:
                sys.stderr.write("%d sequences in alignment\n" % len(thisDict))
            alignDicts.append((thisDict, thisAlign))
            
            endbase = startbase + thisAlign.get_alignment_length() - 1
            thisCharsetString = "charset c%d = %d - %d; [%s]\n" % (num, startbase, endbase, re.sub('.*/', '', options.filenames[num-1]))
            thisCharpartString = "%d:c%d, " % (num, num)
            charsetString += thisCharsetString
            charpartString += thisCharpartString
            startbase = endbase + 1
            num += 1

sys.stderr.write("%d alignments read\n" % len(alignDicts))

#figure out all necessary taxa in the final alignment
allNames = set()
for mydict, alignment in alignDicts:
    allNames |= set(mydict.iterkeys())

sys.stderr.write("%d names across all alignments\n" % len(allNames))

#collect the sequences for each taxon for each alignment, or a dummy string of N's if they are missing from an alignment
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
            finalSeq.id = extract_name_with_regex(options.name_pattern, finalSeq.id)
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

charsetString += charpartString 
charsetString += ";\nend;\n"
sys.stdout.write("%s\n" % charsetString)

