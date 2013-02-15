#!/usr/bin/env python
import sys
import re
import argparse

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from dzutils import ParsedSequenceDescription
from dzutils import read_from_file_or_pickle


'''
#this was just some experimentation with making my own hashes
#it turns out that __cmp__ also needs to be overridden
class my_int(int):
    def __hash__(self):
        val = 42
        print "   (in __hash__: %d)" % val
        return val
    def __cmp__(self, other):
        return False

lint = [ 1, 2, 2 ]
my_lint = [ my_int(i) for i in lint ]

print "lint hash:", lint[0].__hash__()
print "my_lint hash:", my_lint[0].__hash__()

sint = set( lint )
print "about to make set"
my_sint = set (my_lint)

print sint
print my_sint
exit()

#s = set( [ my_int(1), my_int(2), my_int(2) ] )
print s
'''

class my_SeqRecord(SeqIO.SeqRecord, object):
    '''
    def __init__(self, other=None):
        if other is not None:
            self = copy.deepcopy(other)
    '''
    def __hash__(self):
        val = hash(str(self.format('fasta')))
        #print "hashing!: %s" % val
        return val
    def __cmp__(self, other):
        return False

#def add_hash(an_instance):
#    an_instance.__hash__ = new.instancemethod(my_seqrecord_hash, an_instance, an_instance.__class__)

#add_hash(SeqIO.SeqRecord)

#use argparse module to parse commandline input
parser = argparse.ArgumentParser(description='extract sequences from a fasta file')

mutExclGroup = parser.add_mutually_exclusive_group()

#add possible arguments
parser.add_argument('-v', '--invert-match', dest='invertMatch', action='store_true', default=False,
                    help='invert the sense of the match (default false)')

mutExclGroup.add_argument('-s', '--sort', dest='sortOutput', action='store_true', default=False,
                    help='alphanumerically sort the sequences by name (default false)')

mutExclGroup.add_argument('-sc', '--sort-coord', dest='sortOutputByCoord', action='store_true', default=False,
                    help='sort sequences by coordinate if embedded in description (default false)')

parser.add_argument('-p', '--pickle', dest='usePickle', action='store_true', default=False,
                    help='read and write pickled sequence files (<seqfile>.pickle, sometimes faster, default False)')

parser.add_argument('-ms', '--match-sequence', dest='matchSequence', action='store_true', default=False,
                    help='search through the actual sequences for a match, rather than the sequence names')

parser.add_argument('--range', dest='baseRange', nargs=2, type=int, default=[1, -1], metavar=('startbase', 'endbase'),
                    help='range of alignment positions to output, start at 1, last position included, -1 for end')

parser.add_argument('-f', '--patternfile', dest='patternFile', type=str, default=None, 
                    help='file from which to read patterns (you must still pass a pattern on the command line, which is ignored)')

parser.add_argument('pattern',
                    help='a quoted regular expression to search sequence names for')

parser.add_argument('filenames', nargs='*', default=[], 
                    help='a list of filenames to search')

#now process the command line
options = parser.parse_args()

log = sys.stderr

startBase = options.baseRange[0] - 1
endBase = options.baseRange[1]
if options.patternFile is not None:
    log.write('reading patterns from file %s ...\n' % options.patternFile)
    pf = open(options.patternFile, 'rb')
    seqPatterns = [ line.strip() for line in pf ]
    log.write('patterns: %s' % str(seqPatterns))
else:
    seqPatterns = [options.pattern]
seqFiles = options.filenames

compiledPats = []
for pat in seqPatterns:
    try:
        cpat = re.compile(pat)
        compiledPats.append(cpat)
    except StandardError:
        log.write("problem compiling regex pattern %s\n" % pat)
        exit(1)

log.write("Parsing sequence files %s ...\n" % str(seqFiles))

#read the recs from all of the files
recList = []
recSet = set()
for oneSeqFile in seqFiles:
    try:
        #allRecs = [ rec for rec in SeqIO.parse(oneSeqFile, "fasta", alphabet=IUPAC.ambiguous_dna) ]
        #allRecs = set([ rec for rec in SeqIO.parse(oneSeqFile, "fasta", alphabet=IUPAC.ambiguous_dna) ])
        #recList.extend([ rec for rec in SeqIO.parse(oneSeqFile, "fasta", alphabet=IUPAC.ambiguous_dna) ])
        #recList |= set([ rec for rec in SeqIO.parse(oneSeqFile, "fasta", alphabet=IUPAC.ambiguous_dna) ])
        
        if not options.usePickle:
            recList.extend( rec for rec in SeqIO.parse(oneSeqFile, "fasta", alphabet=IUPAC.ambiguous_dna) )
        else:
            recList.extend( rec for rec in read_from_file_or_pickle(oneSeqFile, oneSeqFile + '.pickle', SeqIO.parse, "fasta", alphabet=IUPAC.ambiguous_dna) )
      
        #recSet |= set([ rec for rec in SeqIO.parse(oneSeqFile, "fasta", alphabet=IUPAC.ambiguous_dna) ])
    except IOError:
        log.write("error reading file %s!" % oneSeqFile)
        exit(1)
    
    #if sortOutput:
    #    allRecs.sort(key=lambda x:x.name)

recSet = set( my_SeqRecord(rec.seq, description=rec.description, name=rec.name, id =rec.id) for rec in recList )

if options.invertMatch:
    matchedRecs = set(recSet)
else:
    matchedRecs = set()
for cpat in compiledPats:
    for seq in recSet:
        if options.matchSequence:
            match = cpat.search(seq.seq.data)
        else:
            match = cpat.search(seq.description)
        if match is not None:
            if options.invertMatch:
                if seq in matchedRecs:
                    matchedRecs.remove(seq)
            else:
                matchedRecs.add(seq)
        '''     
        if ( match is not None and not invertMatch ) or ( match is None and invertMatch ):
            #element -1 is the last element, but slicing includes up to but not including the second
            #value.  So, leave it out to get up to the actual end
            if endBase == -1:
                matchedRecs.append(seq[startBase:])
            else:
                matchedRecs.append(seq[startBase:endBase])
        '''

if endBase == -1:
    preparedRecs = [ seq[startBase:] for seq in matchedRecs ]
else:
    preparedRecs = [ seq[startBase:endBase] for seq in matchedRecs ]

#print ParsedSequenceDescription(preparedRecs[0].description)

if options.sortOutput:
    log.write("sorting by sequence name\n")
    preparedRecs.sort(key=lambda x:x.name)
elif options.sortOutputByCoord:
    log.write("sorting by sequence start coordinate\n")
    for rec in preparedRecs:
        rec.coord_start  = ParsedSequenceDescription(rec.description).coord_start
    preparedRecs.sort(key=lambda x:x.coord_start)

if matchedRecs:
    log.write("matched %d sequences in %s\n" % (len(preparedRecs), str(seqFiles)))
    SeqIO.write(preparedRecs, sys.stdout, "fasta")
else:
    log.write("no sequences matched in %s!\n" % str(seqFiles))


