#!/usr/bin/env python

import sys

if __name__ == "__main__":
    import doctest
    doctest.testmod()

class BlinkCluster:
    def __init__(self, num, members=[]):
        self.number = num
        self.cluster_members = members
    def add(self, member):
        self.cluster_members.append(member)
    def __len__(self):
        return len(self.cluster_members)


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
            my_output("problem reading line %s of blink.out\n" % (str(line)), logfile)
            my_output("expecting lines with only:\ncluster# seqname\n", logfile)
            exit(1)
    return allClusters

class HitList:
    def __init__(self, hits):
        #self.hitlist = sets.Set(hits)
        self.hitlist = set(hits)
        self.uniqueNames = None
    
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
            for hit in self.hitlist:
                if not hit[0] in self.uniqueNames:
                    self.uniqueNames[hit[0]] = count
                    count +=1
                if not hit[1] in self.uniqueNames:
                    self.uniqueNames[hit[1]] = count
                    count +=1
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
    

