#!/usr/bin/env python

from re import search, sub
from dzutils import ParsedSequenceDescription
from dzutils import CoordinateSet
from dzutils import extract_core_filename
#from dzutils import *

class Oryza(object):
    def __init__(self):
        #self.taxon_names = [ 'O. sativaj AA', 'O. sativai AA', 'O. barthii AA', 'O. brachyantha FF', 'O. glaberrima AA', 'O. glaberrimaF AA', 'O. glaberrimaM AA', 'O. glumaepatula AA', 'O. meridionalis AA', 'O. minuta BB', 'O. minuta CC', 'O. nivara AA', 'O. officinalis CC', 'O. punctata BB', 'O. rufipogon AA' ]
        #self.short_names = [ 'OsatjAA', 'OsatiAA', 'ObartAA', 'ObracFF', 'OglabAA', 'OglaFAA', 'OglaMAA', 'OglumAA', 'OmeriAA', 'OminuBB', 'OminuCC', 'OnivaAA', 'OoffiCC', 'OpuncBB', 'OrufiAA' ]
        self.taxon_names = [ 'O. barthii AA', 'O. brachyantha FF', 'O. glaberrima AA', 'O. glaberrimaF AA', 'O. glaberrimaM AA', 'O. glumaepatula AA', 'O. meridionalis AA', 'O. minuta BB', 'O. minuta CC', 'O. nivara AA', 'O. officinalis CC', 'O. punctata BB', 'O. rufipogon AA', 'O. sativai AA', 'O. sativaj AA' ]
        self.short_names = [ 'ObartAA', 'ObracFF', 'OglabAA', 'OglaFAA', 'OglaMAA', 'OglumAA', 'OmeriAA', 'OminuBB', 'OminuCC', 'OnivaAA', 'OoffiCC', 'OpuncBB', 'OrufiAA', 'OsatiAA', 'OsatjAA' ]
        self.short_to_long = dict( [ (self.short_names[n], self.taxon_names[n]) for n in range(0, len(self.taxon_names)) ])
        self.long_to_short = dict( [ (self.taxon_names[n], self.short_names[n]) for n in range(0, len(self.taxon_names)) ])
    def short_name_to_long(self, short):
        return self.short_to_long[short]
    def long_name_to_short(self, short):
        return self.long_to_short[short]

oryza = Oryza()


def filter_out_alignments_with_borked_sativa_extractions(toFilter):
    #these are the alignments affected by the sativa mis-extraction problems in /productionOryza2/gramene34_split/alignmentsAndTrees.glabM/ms2006.frac0.5/dag.G1.D2.C4.N5/
    #added here to temporarily easily strip them from various uses of the alignment names
    namesToRemove = [
        '00034.00009.10T.noDupes',
        '00086.00051.11T.noDupes',
        '00146.00032.10T.noDupes',
        '00216.00061.10T.noDupes',
        '00434.00224.11T.noDupes',
        '00477.00041.9T.noDupes',
        '00531.00267.11T.noDupes',
        '00669.00173.10T.noDupes',
        '00675.00175.10T.noDupes',
        '00787.00407.11T.noDupes',
        '01013.00282.10T.noDupes',
        '01043.00474.11T.noDupes',
        '01123.00307.10T.noDupes',
        '01152.00137.9T.noDupes',
        '01261.00059.8T.noDupes',
        '01276.00357.10T.noDupes',
        '01602.00404.10T.noDupes',
        '01639.00100.8T.noDupes',
        '01740.00114.8T.noDupes',
        '01785.00070.7T.noDupes',
        '01835.00060.6T.noDupes'
        ]
    return filter_out_strings_by_pattern(toFilter, namesToRemove)
    

def filter_out_strings_by_pattern(toFilter, patterns):
    if not patterns:
        return toFilter
    filtering = []
    for name in toFilter:
        for patt in patterns:
            if search(patt, name):
                break
        else:
            filtering.append(name)
    return filtering


def rename_sativa_to_oge_standard(name):
    '''standardize the naming format'''
    newName = name.replace('LOC_Os', 'OsatjAA')
    #indica names are like this: BGIOSGA009362
    #the -TA appears in indica gene transcript (mRNA) names. I was stripping this off because it was
    #annoying, but that caused problems because then the gene and mRNA had the same name, and the mRNA 
    #was its own parent.
    if '-TA' in newName:
        newName = newName.replace('BGIOSGA0', 'OsatiAA03gt')
        newName = newName.replace('-TA', '')
    else:
        newName = newName.replace('BGIOSGA0', 'OsatiAA03g')
    #newName = newName.replace('-TA', '')
    #need to make the two glab short names unique
    if 'ORGLA' in name:
        #glaberrima MIPS annotations are named like this: ORGLA03G0400100.1
        #want OglabAA03S_M4001
        newName = newName.replace('ORGLA03G0', 'OglaMAA03S_M')
        newName = sub('00[.](1)', '.\\1', newName)
        newName = sub('00$', '', newName)
    elif 'OglabAA' in name:
        newName = newName.replace('OglabAA', 'OglaFAA')
    return newName.replace('BGIOSIFCE', 'OsatiAA03.')


def extract_all_information_for_seqs_in_alignments(filenames, returnAs='list'):
    '''This script is extracting information from something like the following that I write to the end of the nexus alignments, and returning
    a list of tuples (one per alignment file) with (corefilename, [(seqname, ParsedSequenceDescription)], CoordinateSet)
    Alternatively, if returnAs is 'dict', then return a dict with corefilename keys and (dict(seqname: ParsedSequenceDescription), CoordinateSet) values

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
    if returnAs == 'dict':
        alignments = {}
    else:
        alignments = []

    if isinstance(filenames, str):
        filenames = [ filenames ]
    #work through the files
    for filename in filenames:
        #open a nexus alignment
        alfile = open(filename, 'rb')
        #get the part of the alignment filename that will be identical to part of the treefile name, according to my convention
        coreFilename = extract_core_filename(filename)
        
        #read lines at end of nexus file that give information on sequences, including coordinate
        seqLines = [ line for line in alfile if line.startswith('len') and not 'LOC' in line ]
        #seqLines = [ line for line in file if ((line[0] == 'O' or line[0:3] == 'len') and not 'LOC' in line )]
        try:
            seqDescs = []
            for desc in seqLines:
                found = search('.*(O.*) = (.*)', desc)
                #pull out tuples for normalized taxon names and longer more informative OGE description strings
                seqDescs.append((found.group(1).strip(), found.group(2).strip()))
        except:
            raise RuntimeError('problem parsing file %s' % filename)
        #make a CoordinateSet structure for this alignment file
        coords = CoordinateSet(oryza.taxon_names)
        #parse the description part into my ParsedSequenceDescription data structure
        if returnAs == 'dict':
            parsed = dict([(seq[0], ParsedSequenceDescription(seq[1])) for seq in seqDescs])
            for key, val in parsed.items():
                coords.set_coordinate(key, val.coord_start)
        else:
            parsed = [ (seq[0], ParsedSequenceDescription(seq[1])) for seq in seqDescs ]
            for p in parsed:
                coords.set_coordinate(p[0], p[1].coord_start)
        coords.set_filename(filename)
        #collect a tuple for this alignment with the filename, parsed seq descriptions, and CoordinateSet
        if returnAs == 'dict':
            alignments[str(coreFilename)] = (parsed, coords)
        else:
            alignments.append( (str(coreFilename), parsed, coords) )
    return alignments

