#!/usr/bin/env python

import sys
import re
from dzutils import *

class Oryza:
    def __init__(self):
        self.taxon_names = [ 'O. sativaj AA', 'O. sativai AA', 'O. barthii AA', 'O. brachyantha FF', 'O. glaberrima AA', 'O. glaberrimaF AA', 'O. glaberrimaM AA', 'O. glumaepatula AA', 'O. meridionalis AA', 'O. minuta BB', 'O. minuta CC', 'O. nivara AA', 'O. officinalis CC', 'O. punctata BB', 'O. rufipogon AA' ]
        self.short_names = [ 'OsatjAA', 'OsatiAA', 'ObartAA', 'ObracFF', 'OglabAA', 'OglaFAA', 'OglaMAA', 'OglumAA', 'OmeriAA', 'OminuBB', 'OminuCC', 'OnivaAA', 'OoffiCC', 'OpuncBB', 'OrufiAA' ]
        self.short_to_long = dict( [ (self.short_names[n], self.taxon_names[n]) for n in range(0, len(self.taxon_names)) ])
        self.long_to_short = dict( [ (self.taxon_names[n], self.short_names[n]) for n in range(0, len(self.taxon_names)) ])
    def short_name_to_long(self, short):
        return self.short_to_long[short]
    def long_name_to_short(self, short):
        return self.long_to_short[short]

oryza = Oryza()

def parse_feature_name(feature, errorIsFatal=True):
    '''Return what I'm treating as the name of the SeqFeature, which is 
    stored as one of the arbitrarily named qualifiers, named differently
    for IRGSP and OGE gff's
    '''
    #print feature.qualifiers
    if 'Alias' in feature.qualifiers:
        return feature.qualifiers['Alias'][0]
    elif 'Name' in feature.qualifiers:
        return feature.qualifiers['Name'][0]
    else:
        print 'unable to parse a name!:'
        print feature
        if errorIsFatal:
            exit()

def rename_sativa_to_oge_standard(name):
    '''standardize the naming format'''
    newName = name.replace('LOC_Os', 'OsatjAA')
    #indica names are like this: BGIOSGA009362
    newName = newName.replace('BGIOSGA0', 'OsatiAA03g')
    if 'ORGLA' in name:
        #glaberrima MIPS annotations are named like this: ORGLA03G0400100.1
        #want OglabAA03S_M4001
        newName = newName.replace('ORGLA03G0', 'OglabAA03S_M')
        newName = re.sub('00[.](1)', '.\\1', newName)
        newName = re.sub('00$', '', newName)
    return newName.replace('BGIOSIFCE', 'OsatiAA03.')


def extract_all_information_for_seqs_in_alignments(filenames):
    '''This script is extracting information from something like the following that I write to the end of the nexus alignments, and returning
    a list of tuples (one per alignment file) with (corefilename, [(seqname, ParsedSequenceDescription)]

    [Cluster 82: 10 seq
    len 927    O. sativa AA         = LOC_Os03g29730.1|13103.m03407|CDS expressed protein
    len 954    O. barthii AA        = ObartAA03S_FGT1851 seq=cds; coord=barthii_3s:15315933..15318917:-1; parent_gene=ObartAA03S_FG1851
    len 894    O. brachyantha FF    = ObracFF03S_FGT1583 seq=cds; coord=brachyantha_3s:13893502..13896694:-1; parent_gene=ObracFF03S_FG1583
    len 984    O. glaberrima AA     = OglabAA03S_FGT1853 seq=cds; coord=glaberrima_3s:15719814..15722941:-1; parent_gene=OglabAA03S_FG1853
    len 957    O. minuta BB         = OminuBB03S_FGT1701 seq=cds; coord=minuta_BB_3s:18108842..18112421:-1; parent_gene=OminuBB03S_FG1701
    len 942    O. minuta CC         = OminuCC03S_FGT2016 seq=cds; coord=minuta_CC_3s:22624066..22627315:-1; parent_gene=OminuCC03S_FG2016
    len 897    O. nivara AA         = OnivaAA03S_FGT1699 seq=cds; coord=nivara_3s:14995276..14998560:-1; parent_gene=OnivaAA03S_FG1699
    len 942    O. officinalis CC    = OoffiCC03S_FGT2000 seq=cds; coord=officinalis_3s:23124530..23127781:1; parent_gene=OoffiCC03S_FG2000
    len 999    O. punctata BB       = OpuncBB03S_FGT1941 seq=cds; coord=punctata_3s:19525038..19527936:-1; parent_gene=OpuncBB03S_FG1941
    len 897    O. rufipogon AA      = OrufiAA03S_FGT1588 seq=cds; coord=rufipogon_3s:15199116..15202402:-1; parent_gene=OrufiAA03S_FG1588
    alignment length 1095
    longest sequence 999
    ratio is 0.912329

    the "len ###' bit was added later, and won't appear in earlier alignments
    sativa was later standardized to look like this
    len 3340   O. sativaj AA        = OsatjAA03g29730 seq=gene; coord=Chr3:16936454..16939793:-1
    ]
    '''
    alignments = []

    if isinstance(filenames, str):
        filenames = [ filenames ]
    #work through the files
    for filename in filenames:
        #open a nexus alignment
        file = open(filename, 'rb')
        #get the part of the alignment filename that will be identical to part of the treefile name, according to my convention
        coreFilename = extract_core_filename(filename)
        
        #read lines at end of nexus file that give information on sequences, including coordinate
        seqLines = [ line for line in file if ((line[0] == 'O' or line[0:3] == 'len') and not 'LOC' in line)]
        try:
            seqDescs = []
            for desc in seqLines:
                search = re.search('.*(O.*) = (.*)', desc)
                #pull out tuples for normalized taxon names and longer more informative OGE description strings
                seqDescs.append((search.group(1).strip(), search.group(2).strip()))
        except:
            sys.stderr.write('problem parsing file %s' % filename)
            exit()
        #parse the description part into my ParsedSequenceDescription data structure
        parsed = [ (seq[0], ParsedSequenceDescription(seq[1])) for seq in seqDescs ]
        #make a CoordinateSet structure for this alignment file
        coords = CoordinateSet(oryza.taxon_names)
        for p in parsed:
            coords.set_coordinate(p[0], p[1].coord_start)
        coords.set_filename(filename)
        #collect a tuple for this alignment with the filename, parsed seq descriptions, and CoordinateSet
        alignments.append( (str(coreFilename), parsed, coords) )
    return alignments

