#!/usr/bin/env python

import sys
import re
import argparse
import copy
import string

from BCBio import GFF
from BCBio.GFF import GFFExaminer
from Bio.SeqRecord import SeqRecord

def sort_feature_list(recList):
    '''allow the passed object to be either a list of SeqRecords, which will be sorted based on their
    first feature name, or a single SeqRecord, with a list of features that is to be sorted.
    Prefer the qualifier 'Alias', which is in IRGSP and properly maintains ordering there
    '''
    if isinstance(recList, list):
        qual = 'ID'
        if 'Alias' in recList[0].features[0].qualifiers:
            qual = 'Alias'
        try:
            recList.sort(key=lambda rec:rec.features[0].qualifiers[qual])
        except KeyError:
            sys.stderr.write('ERROR qualifier %s not found\n' % qual)
            for feat in recList.features:
                if qual not in feat.qualifiers:
                    sys.stderr.write('feature:\n%s' % feat)
            exit(1)

   elif isinstance(recList, SeqRecord):
        qual = 'ID'
        if 'Alias' in recList.features[0].qualifiers:
            qual = 'Alias'
        try:
            recList.features.sort(key=lambda feat:feat.qualifiers[qual])
        except KeyError:
            sys.stderr.write('ERROR qualifier %s not found\n' % qual)
            for feat in recList.features:
                if qual not in feat.qualifiers:
                    sys.stderr.write('feature:\n%s' % feat)
            exit(1)

parser = argparse.ArgumentParser(description='extract records from a gff file')

#add possible arguments
parser.add_argument('-v', '--invert-match', dest='invertMatch', action='store_true', default=False,
                    help='invert the sense of the match (default false)')

parser.add_argument('-s', '--sort', dest='sortOutput', action='store_true', default=False,
                    help='alphanumerically sort the sequences by name (default false)')

parser.add_argument('-f', '--patternfile', dest='patternFile', type=str, default=None, 
                    help='file from which to read patterns (you must still pass a pattern on the command line, which is ignored)')

parser.add_argument('pattern',
                    help='a quoted regular expression to search sequence names for')

parser.add_argument('filenames', nargs='*', default=[], 
                    help='a list of filenames to search (none for stdin)')

#now process the command line
parsed = parser.parse_args()

invertMatch = parsed.invertMatch
if parsed.patternFile is not None:
    sys.stderr.write('reading patterns from file %s ...\n' % parsed.patternFile)
    pf = open(parsed.patternFile, 'rb')
    seqPatterns = [ line.strip() for line in pf ]
    sys.stderr.write('patterns: %s\n' % str(seqPatterns))
else:
    seqPatterns = [parsed.pattern]
gffFiles = parsed.filenames
sortOutput = parsed.sortOutput

compiledPats = []
for pat in seqPatterns:
    try:
        cpat = re.compile(pat)
        compiledPats.append(cpat)
    except:
        sys.stderr.write("problem compiling regex pattern %s\n" % pat)
        exit(1)

sys.stderr.write("Parsing gff files %s ...\n" % str(gffFiles))

#this is only dealing with a single gff file at the moment
in_handle = open(gffFiles[0])

allNewRecs = []
#get all of the features for the records (in this case 1, the whole chrom), then we'll get what we want below
#for rec in GFF.parse(in_handle, limit_info=limit_info, base_dict=seq_dict):
for rec in GFF.parse(in_handle):
    sys.stderr.write("%s : %d toplevel features \n" % ( rec.name,  len(rec.features)))
    num = 1
    
    #loop over the features - these will usually be genes
    startFeat = 0
    if len(rec.features) > 0:
        #create a new top level record that can then be repopulated with features and written to a new gff
        newRec = copy.deepcopy(rec)
        newRec.features = []
        
        #if rec is something above the gene/mRNA level, start downstream of it
        if string.lower(rec.features[0].type) in [ 'chromosome', 'contig', 'scaffold' ]:
            startFeat = 1

        if invertMatch:
            matchedFeats = set(rec.features[startFeat:])
        else:
            matchedFeats = set()
        for cpat in compiledPats:
            #loop over features of the rec
            for feat in rec.features[startFeat:]:
                hit = False
                for qual in feat.qualifiers.items():
                    #for a given qualifier the value is a list of strings
                    for it in qual[1]:
                        match = cpat.search(it)
                        if match is not None:
                            hit = True
                            break
                    if hit:
                        if invertMatch:
                            if feat in matchedFeats:
                                matchedFeats.remove(feat)
                        else:
                            matchedFeats.add(feat)
                        break
        newRec.features.extend(list(matchedFeats))
        sort_feature_list(newRec)
    allNewRecs.append(newRec)
         
if sortOutput:
    sort_feature_list(allNewRecs)

#gffOutFilename = re.sub('^.*\/', '', parsed.gffFilename) + '.extracted'
#gffOut = open(gffOutFilename, "w")

gffOut = sys.stdout
GFF.write( allNewRecs , gffOut)

