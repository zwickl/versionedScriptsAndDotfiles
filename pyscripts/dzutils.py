#!/usr/bin/env python

import sys
import re

if __name__ == "__main__":
    import doctest
    doctest.testmod()

class BlinkCluster:
    def __init__(self, num, members=[], mapping=None):
        self.number = num
        if mapping is None:
            self.cluster_members = members
        else:
            self.cluster_members = []
            for member in members:
                self.cluster_members.append(mapping[member])

    def add(self, member):
        self.cluster_members.append(member)
    def __len__(self):
        return len(self.cluster_members)
    def output(self, stream=sys.stdout):
        for mem in self.cluster_members:
            stream.write('%d\t%s\n' % (self.number, mem))

def parse_mcl_output(filename):
    '''read MCL output, which looks like the below, return a list of BlinkClusters
    output is simply one line per cluster:
    OminuCC03S_FGT1789  OminuCC03S_FGT1792  OoffiCC03S_FGT1798  OminuBB03S_FGT0728  OminuBB03S_FGT0727  OoffiCC03S_FGT1799  OminuCC03S_FGT1791  OminuCC03S_FGT1790   
    '''
    mclOut = open(filename, "rb")
    lines = [ l.split() for l in mclOut ]
    allClusters = []

    num = 1
    for cluster in lines:
        try:
            allClusters.append(BlinkCluster(num, cluster))
            num = num + 1
        except:
            exit("problem converting mcl cluster to blink format")
    return allClusters

def parse_blink_output(filename):
    '''read blink output, which looks like the below, return a list of BlinkClusters
    This indicates cluster 0 with one member, cluster 1 with 4, etc.
    0	ObartAA03S_FGT0005
    1	OglabAA03S_FGT0268
    1	OrufiAA03S_FGT0184
    1	OglabAA03S_FGT0269
    1	ObartAA03S_FGT0026
    2	OrufiAA03S_FGT0182
    2	OglabAA03S_FGT0266
    2	OminuCC03S_FGT0238
    '''
    blinkOut = open(filename, "rb")
    lines = [ l.split() for l in blinkOut ]
    allClusters = []
    thisCluster = []
    curNum = 0
    first = True
    for line in lines:
        try:
            if len(line) != 0 and len(line) != 2:
                raise Exception
            if len(line) == 0:
                thisLineNum = -1
            else:
                thisLineNum = int(line[0])
                #if the number is a float
                if str(thisLineNum) != line[0]:
                    raise Exception
            if first is True:
                curNum = thisLineNum
                first = False
                
            #if we haven't hit a blank line (presumably the end of the file)
            #and this line has the same cluster number as previous ones, append it
            if thisLineNum == curNum:
                thisCluster.append(line[1])
            #if we did hit a blank line, or this is the last line in the file
            #or we found a new cluster number, store the current cluster
            if line == lines[-1] or thisLineNum != curNum:
                toAppend = thisCluster
                #allClusters.append([curNum, toAppend])
                allClusters.append(BlinkCluster(curNum, toAppend))
                if thisLineNum >= 0 and int(line[0]) != curNum:
                    curNum = thisLineNum
                    thisCluster = []
                    thisCluster.append(line[1])
        except:
            print "problem reading line %s of blink.out\n" % (str(line))
            print "expecting lines with only:\ncluster# seqname\n"
            #my_output("problem reading line %s of blink.out\n" % (str(line)), logfile)
            #my_output("expecting lines with only:\ncluster# seqname\n", logfile)
            exit(1)
    return allClusters

'''
def blink_cluster_from_clique(thisClust, maxClique, mapping=None):
    if mapping is not None:
        new_members = []
        for member in clique:
            if mapping is not None:
                new_members.append(mapping[member]
            else:
                new_members.append(member)
'''

class HitList:
    def __init__(self, hits):
        #self.hitlist = sets.Set(hits)
        self.hitlist = set(hits)
        self.uniqueNames = None
        self.numbersToNames = None
    
    def __len__(self):
        return len(self.hitlist)

    def get_sublist_by_query_names(self, names):
        #subset = sets.Set()
        subset = set()
        for hit in self.hitlist:
            for name in names:
                if name in hit[0]:
                    subset.add(hit)
        #return subset
        return HitList(subset)
        
    def get_sublist_by_hit_names(self, names):
        subset = set()
        for hit in self.hitlist:
            for name in names:
                if name in hit[1]:
                    subset.add(hit)
        return HitList(subset)

    def get_sublist_by_query_or_hit_names(self, names):
        subset = self.get_sublist_by_query_names(names)
        return subset.union(self.get_sublist_by_hit_names(names))
   
    def union(self, others):
        return HitList(self.hitlist.union(others.hitlist))

    def unique_names(self):
        if self.uniqueNames is None:
            count = 1
            self.uniqueNames = {}
            self.numbersToNames = {}
            for hit in self.hitlist:
                for name in hit:
                    if not name in self.uniqueNames:
                        self.uniqueNames[name] = count
                        self.numbersToNames[count] = name
                        count +=1
            if len(self.uniqueNames) != len(self.numbersToNames):
                exit("problem mapping hit names to numbers")
        return self.uniqueNames.keys()
        '''
            self.uniqueNames = set()
            for hit in self.hitlist:
                self.uniqueNames.add(hit[0])
                self.uniqueNames.add(hit[1])
        return self.uniqueNames
        '''

    def output(self, stream=sys.stdout):
        for hit in self.hitlist:
            stream.write("%s\t%s\n" % hit)

    def get_list_numbers_for_names(self):
        if self.uniqueNames is None:
            self.unique_names()
        numSet = set()
        for hit in self.hitlist:
            set.add((self.uniqueNames[hit[0]], self.uniqueNames[hit[1]]))

    def output_for_dfmax(self, stream=sys.stdout):
        if self.uniqueNames is None:
            self.unique_names()

        stream.write('p EDGE %s %s\n' % (str(len(self.uniqueNames)), str(len(self.hitlist))))
        for hit in self.hitlist:
            stream.write('e %s %s\n' % (self.uniqueNames[hit[0]], self.uniqueNames[hit[1]]))

def parse_hits_file(filename):
    hitsFile = open(filename, "rb")
    lines = [ tuple(l.split()) for l in hitsFile ]
    #return HitList(lines).get_sublist_by_query_names(['LOC'])
    return HitList(lines)
   
def make_dictionary_from_gff_arbitrary_field(string):
    dict = {}
    fields = string.split(';')
    for f in fields:
        sep = f.split('=')
        try:
            dict[sep[0]] = sep[1]
        except:
            exit("problem reading field, %s" % sep)
    return dict

class ParsedSequenceDescription:
    def __init__(self, description=None, gff=None):
        '''
        OGE gff lines:
        >>> print ParsedSequenceDescription(gff='barthii_3s  ensembl gene    735328  738974  .   -   .   ID=ObartAA03S_FG0284;Name=ObartAA03S_FG0284;biotype=protein_coding').name
        ObartAA03S_FG0284
        >>> print ParsedSequenceDescription(gff='barthii_3s  ensembl mRNA    735328  738974  .   -   .   ID=ObartAA03S_FGT0284;Parent=ObartAA03S_FG0284;Name=ObartAA03S_FGT0284;biotype=protein_coding').name
        ObartAA03S_FGT0284
        >>> print ParsedSequenceDescription(gff='barthii_3s  ensembl CDS 738588  738974  .   -   0   Parent=ObartAA03S_FGT0284;Name=CDS.1419').name
        ObartAA03S_FGT0284
        
        OGE fasta sequence names
        >>> print ParsedSequenceDescription(description='ObartAA03S_FGT0284 seq=cds; coord=barthii_3s:735328..738974:-1; parent_gene=ObartAA03S_FG0284').name
        ObartAA03S_FGT0284
        >>> print ParsedSequenceDescription(description='>ObartAA03S_FG0284 seq=gene; coord=barthii_3s:735328..738974:-1').name
        ObartAA03S_FG0284
        
        IRGSP gff
        >>> print ParsedSequenceDescription(gff='Chr3    MSU_osa1r6  gene    934041  938212  .   -   .   ID=13103.t00151;Name=proteasome%20subunit%2C%20putative%2C%20expressed;Alias=LOC_Os03g02540').name
        LOC_Os03g02540
        >>> print ParsedSequenceDescription(gff='Chr3    MSU_osa1r6  mRNA    934041  938212  .   -   .   ID=13103.m00215;Parent=13103.t00151;Alias=LOC_Os03g02540.1').name
        LOC_Os03g02540.1
        >>> print ParsedSequenceDescription(gff='Chr3    MSU_osa1r6  CDS 937701  938087  .   -   0   Parent=13103.m00215').name
        13103.m00215
        
        IRGSP fasta sequence names
        >>> print ParsedSequenceDescription(description='>LOC_Os03g02540.1|13103.m00215|CDS proteasome subunit, putative, expressed').name
        LOC_Os03g02540.1
        >>> print ParsedSequenceDescription(description='>LOC_Os03g02540|13103.t00151|unspliced-genomic proteasome subunit, putative, expressed').name
        LOC_Os03g02540
        
        '''

        if gff is not None and description is not None:
            exit("pass either a gff or description string, not both")

        self.name = None
        self.ID = None
        self.description = None
        self.type = None
        self.molecule = None
        self.coord_start = None
        self.coord_end = None
        self.strand = None
        self.frame = None
        self.parent = None
        '''
        #############################
        OGE:

        CDS fasta description:
        >ObartAA03S_FGT0284 seq=cds; coord=barthii_3s:735328..738974:-1; parent_gene=ObartAA03S_FG0284

        gene fasta description
        >ObartAA03S_FG0284 seq=gene; coord=barthii_3s:735328..738974:-1

        gff for this looks like:
        barthii_3s  ensembl gene    735328  738974  .   -   .   ID=ObartAA03S_FG0284;Name=ObartAA03S_FG0284;biotype=protein_coding
        barthii_3s  ensembl mRNA    735328  738974  .   -   .   ID=ObartAA03S_FGT0284;Parent=ObartAA03S_FG0284;Name=ObartAA03S_FGT0284;biotype=protein_coding
        barthii_3s  ensembl intron  737785  738587  .   -   .   Parent=ObartAA03S_FGT0284;Name=intron.1408
        barthii_3s  ensembl intron  736782  737463  .   -   .   Parent=ObartAA03S_FGT0284;Name=intron.1409
        barthii_3s  ensembl intron  736578  736664  .   -   .   Parent=ObartAA03S_FGT0284;Name=intron.1410
        barthii_3s  ensembl intron  735818  736491  .   -   .   Parent=ObartAA03S_FGT0284;Name=intron.1411
        barthii_3s  ensembl intron  735508  735582  .   -   .   Parent=ObartAA03S_FGT0284;Name=intron.1412
        barthii_3s  ensembl exon    738588  738974  .   -   .   Parent=ObartAA03S_FGT0284;Name=ObartAA03S_FG0284.exon1
        barthii_3s  ensembl exon    737464  737784  .   -   .   Parent=ObartAA03S_FGT0284;Name=ObartAA03S_FG0284.exon2
        barthii_3s  ensembl exon    736665  736781  .   -   .   Parent=ObartAA03S_FGT0284;Name=ObartAA03S_FG0284.exon3
        barthii_3s  ensembl exon    736492  736577  .   -   .   Parent=ObartAA03S_FGT0284;Name=ObartAA03S_FG0284.exon4
        barthii_3s  ensembl exon    735583  735817  .   -   .   Parent=ObartAA03S_FGT0284;Name=ObartAA03S_FG0284.exon5
        barthii_3s  ensembl exon    735328  735507  .   -   .   Parent=ObartAA03S_FGT0284;Name=ObartAA03S_FG0284.exon6
        barthii_3s  ensembl CDS 738588  738974  .   -   0   Parent=ObartAA03S_FGT0284;Name=CDS.1419
        barthii_3s  ensembl CDS 737464  737784  .   -   0   Parent=ObartAA03S_FGT0284;Name=CDS.1420
        barthii_3s  ensembl CDS 736665  736781  .   -   0   Parent=ObartAA03S_FGT0284;Name=CDS.1421
        barthii_3s  ensembl CDS 736492  736577  .   -   0   Parent=ObartAA03S_FGT0284;Name=CDS.1422
        barthii_3s  ensembl CDS 735583  735817  .   -   2   Parent=ObartAA03S_FGT0284;Name=CDS.1423
        barthii_3s  ensembl CDS 735328  735507  .   -   0   Parent=ObartAA03S_FGT0284;Name=CDS.1424
        ######################################
        IRGSP:
        CDS fasta description:
        >LOC_Os03g02540.1|13103.m00215|CDS proteasome subunit, putative, expressed
        
        gene fasta description:
        >LOC_Os03g02540|13103.t00151|unspliced-genomic proteasome subunit, putative, expressed

        gff for this looks like:
        Chr3    MSU_osa1r6  gene    934041  938212  .   -   .   ID=13103.t00151;Name=proteasome%20subunit%2C%20putative%2C%20expressed;Alias=LOC_Os03g02540
        Chr3    MSU_osa1r6  mRNA    934041  938212  .   -   .   ID=13103.m00215;Parent=13103.t00151;Alias=LOC_Os03g02540.1
        Chr3    MSU_osa1r6  five_prime_UTR  938088  938212  .   -   .   Parent=13103.m00215
        Chr3    MSU_osa1r6  CDS 937701  938087  .   -   0   Parent=13103.m00215
        Chr3    MSU_osa1r6  CDS 936586  936906  .   -   0   Parent=13103.m00215
        Chr3    MSU_osa1r6  CDS 935789  935905  .   -   0   Parent=13103.m00215
        Chr3    MSU_osa1r6  CDS 935616  935701  .   -   0   Parent=13103.m00215
        Chr3    MSU_osa1r6  CDS 934712  934946  .   -   2   Parent=13103.m00215
        Chr3    MSU_osa1r6  CDS 934457  934636  .   -   0   Parent=13103.m00215
        Chr3    MSU_osa1r6  three_prime_UTR 934041  934456  .   -   .   Parent=13103.m00215

        '''
        #this gff reading attempt was aborted, I think
        if gff is not None:
            #gff's are generally the same, except for random junk in the last field.  
            #g = gff.replace(";", " ")
            split_gff = gff.split()
            #So, we end up with
            # 0          1        2       3        4     5   6   7      8  
            #Chr3    MSU_osa1r6  gene    3465    5944    .   +   .   ID=13103.t05666;Name=expressed%20protein;Alias=LOC_Os03g01008
            
            # 0             1     2     3       4    5   6   7      8                    
            #nivara_3s   ensembl CDS 1001307 1001549 .   -   0   Parent=OnivaAA03S_FGT0137;Name=CDS.26236
            self.molecule = split_gff[0]
            self.program = split_gff[1]
            self.type = split_gff[2]
            self.coord_start = split_gff[3]
            self.coord_end = split_gff[4]
            self.score = split_gff[5]
            self.strand = split_gff[6]

            self.various = make_dictionary_from_gff_arbitrary_field(split_gff[8])

            if self.type in [ 'gene', 'mRNA' ] :
                if 'Alias' in self.various:
                    self.name = self.various['Alias']
                elif 'Name' in self.various:
                    self.name = self.various['Name']
            else:
                try:
                    self.name = self.various['Parent']
                except:
                    exit("problem extracting name from %s" % self.various)
        else:
            desc = description.replace(";", " ")
            #desc = desc.replace(":", " ")
            split_desc = desc.split()
            
            #get the name
            self.name = split_desc[0].replace(">", "")

            #get the sequence type
            if not "seq" in split_desc[1]:
                exit("no \"seq\" found in %s" % description)
            match = re.search("seq=(.*)", split_desc[1])
            if match is None:
                exit("problem matching \"seq\" in %s" % description)
            self.type = match.group(1)
            
            #get molecule, coordinates and strand
            if not "coord" in split_desc[2]:
                exit("no \"coord\" found in %s" % description)
            match = re.search("coord=(.*):(.*):(.*)", split_desc[2])
            if match is None:
                exit("problem matching \"coord\" in %s" % description)
            self.molecule = match.group(1)
            coordSplit = match.group(2).split(".")
            (self.coord_start, self.coord_end) = (coordSplit[0], coordSplit[2])
            self.strand = match.group(3)

            #get parent gene
            if len(split_desc) > 3:
#made some hasty changes here - check it through at some point
                if not "parent_gene" in split_desc[3]:
                    self.parent = 'none'
                    #exit("no \"parent_gene\" found in %s" % description)
                else:
                    match = re.search("parent_gene=(.*)", split_desc[3])
                    if match is None:
                        exit("problem matching \"parent_gene\" in %s" % description)
                    self.parent = match.group(1)
            else:
                #if this is a full gene, it has no parent
                self.parent = 'none'

    def output(self):
        print "name", self.name
        print "type", self.type
        print "molecule", self.molecule
        print "start", self.coord_start
        print "end", self.coord_end
        print "strand", self.strand
        print "parent", self.parent

if __name__ == "__main__":
    import doctest
    doctest.testmod()

