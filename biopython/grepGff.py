#!/usr/bin/env python
import sys
import re
import argparse
import copy
import string

from BCBio import GFF
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from dzutils import read_from_file_or_pickle
from dzbiopython import  sort_feature_list, sort_feature_list_by_coordinate

parser = argparse.ArgumentParser(description='extract records from a gff file')

mutExclGroup = parser.add_mutually_exclusive_group()

#add possible arguments
parser.add_argument('-v', '--invert-match', action='store_true', default=False,
                    help='invert the sense of the match (default false)')

mutExclGroup.add_argument('-s', '--sort', action='store_true', default=False,
                    help='alphanumerically sort the sequences by name (default false)')

mutExclGroup.add_argument('-sc', '--sort-coord', action='store_true', default=False,
                    help='sort the sequences by start coordinatre (default false)')

parser.add_argument('-p', '--pickle', action='store_true', default=False,
                    help='read and write pickled sequence files (<gfffile>.pickle, sometimes faster, default False)')

parser.add_argument('-f', '--patternfile', type=str, default=None, 
                    help='file from which to read patterns (you must still pass a pattern on the command line, which is ignored)')

parser.add_argument('pattern',
                    help='a quoted regular expression to search sequence names for')

parser.add_argument('--range', nargs=2, type=int, default=[1, -1], metavar=('startbase', 'endbase'),
                    help='only output annotations entirely within these coordinates, start at 1, last position included, -1 for end')

parser.add_argument('filenames', nargs='*', default=[], 
                    help='a list of filenames to search (none for stdin)')

#now process the command line
options = parser.parse_args()

log = sys.stderr

start_base = options.range[0] - 1
end_base = options.range[1]

if options.patternfile:
    log.write('reading patterns from file %s ...\n' % options.patternfile)
    seqPatterns = [ line.strip() for line in open(options.patternfile, 'rb') ]
    log.write('patterns: %s\n' % str(seqPatterns))
else:
    seqPatterns = [options.pattern]

compiledPats = []
for pat in seqPatterns:
    try:
        compiledPats.append(re.compile(pat))
    except StandardError:
        log.write("problem compiling regex pattern %s\n" % pat)
        exit(1)

log.write("Parsing gff files %s ...\n" % str(options.filenames))

#this is only dealing with a single gff file at the moment
#in_handle = open(options.filenames[0])

allNewRecs = []
#get all of the features for the records (in this case 1, the whole chrom), then we'll get what we want below
#for rec in GFF.parse(in_handle, limit_info=limit_info, base_dict=seq_dict):
#for rec in GFF.parse(in_handle):
#for rec in GFF.parse(options.filenames[0]):
if not options.pickle:
    gffRecs = GFF.parse(options.filenames[0])
else:
    gffRecs = read_from_file_or_pickle(options.filenames[0], options.filenames[0] + '.pickle', GFF.parse)

for rec in gffRecs:
    log.write("%s : %d toplevel features \n" % ( rec.name,  len(rec.features)))
    num = 1
    
    #loop over the features - these will usually be genes
    startFeat = 0
    if rec.features:
        #create a new top level record that can then be repopulated with features and written to a new gff
        newRec = copy.deepcopy(rec)
        newRec.features = []
        newRec.annotations = {}
        
        #if rec is something above the gene/mRNA level, start downstream of it
        if string.lower(rec.features[0].type) in [ 'chromosome', 'contig', 'scaffold' ]:
            startFeat = 1

        matchedFeats = set(rec.features[startFeat:]) if options.invert_match else set()
        for cpat in compiledPats:
            #loop over features of the rec
            for feat in rec.features[startFeat:]:
                hit = False
                for qual in feat.qualifiers.items():
                    #for a given qualifier the value is a list of strings
                    for it in qual[1]:
                        match = cpat.search(it)
                        if match:
                            hit = True
                            break
                    if hit:
                        if options.invert_match:
                            #with invert_match the feat is removed from the set
                            if feat in matchedFeats:
                                matchedFeats.remove(feat)
                        else:
                            #otherwise it is added
                            matchedFeats.add(feat)
                        #either way we can stop looping over patterns
                        break
        #whatever is in matchedFeats (possibly nothing) can now be added to 
        #the new rec if it is in the required base range
        rangeFilteredFeats = set(matchedFeats)
        for feat in matchedFeats:
            assert feat.location.start < feat.location.end
            if feat.location.start.position < start_base or ( feat.location.end.position > end_base and end_base > -1 ):
                rangeFilteredFeats.remove(feat)
                log.write('%s matched but not in base range\n' % feat.qualifiers['ID'])
        newRec.features.extend(list(rangeFilteredFeats))
        sort_feature_list(newRec)
    if newRec.features:
        allNewRecs.append(newRec)

if options.sort:
    sort_feature_list(allNewRecs)
elif options.sort_coord:
    sort_feature_list_by_coordinate(allNewRecs)

#gffOutFilename = re.sub('^.*\/', '', options.gffFilename) + '.extracted'
#gffOut = open(gffOutFilename, "w")

gffOut = sys.stdout
GFF.write(allNewRecs, gffOut)

