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
    
class ParsedSequenceDescription:
    def __init__(self, description):
        '''
        this will work on our files with sequence names like this:
        ObartAA03S_FGT0284 seq=cds; coord=barthii_3s:735328..738974:-1; parent_gene=ObartAA03S_FG0284
        but not the sativa genome names, like this
        LOC_Os03g02540.1|13103.m00215|CDS proteasome subunit, putative, expressed
        '''
        if 'LOC' in description:
            print "can't parse descriptions like\n\t%s" % description
            print "only like\n\tObartAA03S_FGT0284 seq=cds; coord=barthii_3s:735328..738974:-1; parent_gene=ObartAA03S_FG0284"
            exit(1)

        desc = description.replace(";", " ")
        #desc = desc.replace(":", " ")
        split_desc = desc.split()
        
        #get the name
        self.name = split_desc[0]

        #get the sequence type
        if not "seq" in split_desc[1]:
            exit("no \"seq\" found in %s" % description)
        match = re.search("seq=(.*)", split_desc[1])
        if match is None:
            exit("problem matching \"seq\" in %s" % description)
        self.type = match.group(1)
        
        #get molecular, coordinates and strand
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
        if not "parent_gene" in split_desc[3]:
            exit("no \"parent_gene\" found in %s" % description)
        match = re.search("parent_gene=(.*)", split_desc[3])
        if match is None:
            exit("problem matching \"parent_gene\" in %s" % description)
        self.parent = match.group(1)

    def output(self):
        print "name", self.name
        print "type", self.type
        print "molecule", self.molecule
        print "start", self.coord_start
        print "end", self.coord_end
        print "strand", self.strand
        print "parent", self.parent


