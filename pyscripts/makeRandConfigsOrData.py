#!/usr/bin/env python
import os
import sys
import string
import StringIO
import random

random.seed()

#EDIT THESE BEFORE RUNNING
NTAX = 11
NCHAR = 2178
EXCL_MIN = 0.5
EXCL_MAX = 0.99
DEL_MIN = 0.0
DEL_MAX = 0.6

'''
THIS IS HOW THEY WERE SET FOR RBCL
NTAX = 500
NCHAR = 1395
EXCL_MIN = 0.5
EXCL_MAX = 0.99
DEL_MIN = 0.8
DEL_MAX = 0.97
'''

if not len(sys.argv) == 3:
    print "usage: <script> [n/c/a/d] (for nuc, codon, aminoacid or data creation) [num reps]"
    sys.exit(1)
elif not (sys.argv[1] == "n" or sys.argv[1] == "c" or sys.argv[1] == "a" or sys.argv[1] == "d"):
    print "usage: <script> [n/c/a/d] (for nuc, codon, aminoacid or data creation) [num reps]"
    sys.exit(1)
else:
    type = sys.argv[1]
    num = int(sys.argv[2])

#CLASS DEFINITIONS
class TextOption:
    def __init__(self, n, s):
        self.name = n
        self.settings = s
    def choose(self):
        return "%s = %s\n" % (self.name, random.choice(self.settings))
    def __str__(self):
        return "name:" + self.name + " options " + [arg for arg in self.settings]
    def rep(self):
        return "name:" + self.name + " options " + str([arg for arg in self.settings])
        
class IntOption:
    def __init__(self, n, min, max):
        self.name = n
        self.minval = min
        self.maxval = max
    def choose(self):
        return "%s = %d\n" % (self.name, random.randint(self.minval, self.maxval))
    def __str__(self):
        return "name:" + self.name + " options " + [arg for arg in self.settings]
    def rep(self):
        #return "name:" + self.name + " range " + self.minval + "-" + self.maxval
        return "name: %s range %d-%d" % (self.name, self.minval, self.maxval)
    
class FloatOption:
    def __init__(self, n, min, max):
        self.name = n
        self.minval = min
        self.maxval = max
    def choose(self):
        return "%s = %f\n" % (self.name, random.uniform(self.minval, self.maxval))
    def __str__(self):
        return "name:" + self.name + " options " + [arg for arg in self.settings]
    def rep(self):
        return "name:" + self.name + " options " + str([arg for arg in self.settings])

class PairedOption:
    def __init__(self, opts):
        self.name = "pair"
        self.pairList = opts
        #self.opt1 = one
        #self.opt2 = two
    def choose(self):
        #return [arg.choose() for arg in self.pairList]
        pair = random.choice(self.pairList)
        #print pair
        return pair[0].choose() + pair[1].choose()

class ModelSection:
    def __init__(self):
        self.entries = []
        
    def set_default_nuc(self):
        self.entries = []
        self.entries.append(TextOption("datatype", ["nucleotide"]))
        self.entries.append(TextOption("ratematrix", ["6rate", "2rate", "6rate", "(a b a b a c)"]))
        self.entries.append(TextOption("statefrequencies", ["estimate", "equal", "empirical"]))
        self.entries.append(PairedOption([[TextOption("ratehetmodel", ["none"]), TextOption("numratecats", ["1"])], [TextOption("ratehetmodel", ["gamma"]), IntOption("numratecats", 2, 6)]]))
        self.entries.append(TextOption("invariantsites", ["none", "estimate"]))

    def set_default_codon(self):
        self.entries = []
        self.entries.append(TextOption("datatype", ["codon"]))
        self.entries.append(TextOption("geneticcode", ["standard"]))
        self.entries.append(TextOption("ratematrix", ["6rate", "2rate", "6rate", "(a b a b a c)"]))
        self.entries.append(TextOption("statefrequencies", ["equal", "empirical", "f1x4", "f3x4"]))
        self.entries.append(PairedOption([[TextOption("ratehetmodel", ["none"]), TextOption("numratecats", ["1"])], [TextOption("ratehetmodel", ["nonsynonymous"]), IntOption("numratecats", 2, 6)]]))
        self.entries.append(TextOption("invariantsites", ["none"]))
        
    def set_default_aa(self):
        self.entries = []
        self.entries.append(TextOption("datatype", ["codon-aminoacid"]))
        self.entries.append(TextOption("ratematrix", ["poisson", "jones", "wag", "dayhoff", "mtmam", "mtrev", "estimate"]))
        self.entries.append(TextOption("statefrequencies", ["estimate", "equal", "empirical", "jones", "wag", "dayhoff", "mtmam", "mtrev"]))
        self.entries.append(PairedOption([[TextOption("ratehetmodel", ["none"]), TextOption("numratecats", ["1"])], [TextOption("ratehetmodel", ["gamma"]), IntOption("numratecats", 2, 6)]]))
        self.entries.append(TextOption("invariantsites", ["none", "estimate"]))
        
    def remove(self, nameStr):
        for i in range(0, len(self.entries)):
            if str(self.entries[i].name) == str(nameStr):
                self.entries.pop(i)
                return i
        return len(self.entries)

    def replace(self, option):
        index = self.remove(option.name)
        self.entries.insert(index, option)

    def choose(self):
        return "".join([e.choose() for e in self.entries])
      
class GeneralSection:
    def __init__(self):
        self.entries = []
        
    def set_default(self):
        self.entries = []
        #self.entries.append(TextOption("datafname", ["rana.tiny.nex"]))
        self.entries.append(TextOption("constraintfile", ["none"]))
        self.entries.append(TextOption("streefname", ["stepwise", "random", "FILE"]))
        self.entries.append(IntOption("attachmentspertaxon", 1, 100))
        self.entries.append(TextOption("ofprefix", ["poo"]))
        self.entries.append(TextOption("randseed", ["-1"]))
        self.entries.append(TextOption("availablememory", ["256", "512", "1024", "2048"] ) )
        self.entries.append(TextOption("logevery", ["1", "10", "100", "1000", "10000000"]) )
        self.entries.append(TextOption("saveevery", ["1", "10", "100", "1000", "10000000"]) )
        self.entries.append(TextOption("refinestart", ["0", "1"]) )
        self.entries.append(TextOption("outputeachbettertopology", ["0", "1"]) )
        self.entries.append(TextOption("outputcurrentbesttopology", ["0", "1"]))
        self.entries.append(TextOption("enforcetermconditions", ["1"]))
        self.entries.append(IntOption("genthreshfortopoterm", 1000, 20000) )
        self.entries.append(FloatOption("scorethreshforterm", 0.001, 1.0) )
        self.entries.append(FloatOption("significanttopochange", 0.001, 1.0) )
        self.entries.append(TextOption("outputphyliptree", ["0", "1"]))
        self.entries.append(TextOption("outputmostlyuselessfiles", ["0", "1"]))
        self.entries.append(TextOption("writecheckpoints", ["0", "1"]))
        self.entries.append(TextOption("restart", ["0"]))
        self.entries.append(TextOption("outgroup", ["1", "2", "1 5", "1 6 2", "1-4", "2-4"]) )
        self.entries.append(IntOption("searchreps", 2, 5))
        self.entries.append(TextOption("collapsebranches", ["0", "1"]))
        self.entries.append(TextOption("outputsitelikelihoods", ["0"]))    
    
    def remove(self, nameStr):
        for i in range(0, len(self.entries)):
            if str(self.entries[i].name) == str(nameStr):
                self.entries.pop(i)
                return i
        return len(self.entries)

    def replace(self, option):
        index = self.remove(option.name)
        self.entries.insert(index, option)
            
    def choose(self):
        return "".join([e.choose() for e in self.entries])

class MasterSection:
    def __init__(self):
        self.entries = []
        
    def set_default(self):
        self.entries = []
        self.entries.append(IntOption("nindivs", 3, 6))
        self.entries.append(TextOption("holdover", ["1"]))
        self.entries.append(TextOption("holdoverpenalty", ["0.0"]))
        
        self.entries.append(FloatOption("selectionintensity", 0.001, 1.0) )
        self.entries.append(IntOption("stopgen", 500000, 500000) )
        self.entries.append(IntOption("stoptime", 14400, 14400) )
        
        self.entries.append(TextOption("startoptprec", ["0.1", "0.5", "1.0", "2.0"]) )
        self.entries.append(TextOption("minoptprec", ["0.0001", "0.001", "0.01", "0.05", "0.1"]) )
        self.entries.append(IntOption("numberofprecreductions", 0, 40) )
        self.entries.append(FloatOption("treerejectionthreshold", 0.01, 500.0) )
        self.entries.append(FloatOption("topoweight", 0.01, 5.0) )
        self.entries.append(FloatOption("modweight", 0.01, 5.0) )
        self.entries.append(FloatOption("brlenweight", 0.01, 5.0) )
        self.entries.append(FloatOption("randnniweight", 0.01, 5.0) )
        self.entries.append(FloatOption("randsprweight", 0.01, 5.0) )
        self.entries.append(FloatOption("limsprweight",  0.01, 5.0) )
        self.entries.append(TextOption("intervallength", ["10", "100", "200", "500"]) )
        self.entries.append(IntOption("intervalstostore", 1, 10) )
        
        self.entries.append(IntOption("limsprrange", 1, 10) )
        self.entries.append(IntOption("meanbrlenmuts", 1, 10) )
        self.entries.append(IntOption("gammashapebrlen", 1, 2000) )
        self.entries.append(IntOption("gammashapemodel", 1, 2000) )
        self.entries.append(TextOption("uniqueswapbias", ["0.001", "0.01", "0.1", "1.0", "2.0"]) )
        self.entries.append(TextOption("distanceswapbias", ["0.001", "0.01", "0.1", "1.0", "2.0"]) )
        
        self.entries.append(FloatOption("resampleproportion", 1.0, 1.0) )
        self.entries.append(PairedOption([[IntOption("bootstrapreps", 1, 2), TextOption("inferinternalstateprobs", ["0"])], [IntOption("bootstrapreps", 0, 0), TextOption("inferinternalstateprobs", ["0", "1"])]]))

    def remove(self, nameStr):
        for i in range(0, len(self.entries)):
            if str(self.entries[i].name) == str(nameStr):
                self.entries.pop(i)
                return i
        return len(self.entries)

    def replace(self, option):
        index = self.remove(option.name)
        self.entries.insert(index, option)

    def choose(self):
        return "".join([e.choose() for e in self.entries])

#RANDOM FUNCTIONS
def get_better_random_set_triples(end, prop):
    chars = range(int(end/3))
    a = int(end*prop/3.0)
    s = random.sample(chars, a)
    s.sort()
    B =  " ".join(["%d-%d" %(3*i+1, 3*i+3) for i in s])
    #print "Exset * poo = %s;" % B
    return B

def get_better_random_set(end, prop):
    chars = range(end)
    s = random.sample(chars, int(end * prop))
    s.sort()
    B =  " ".join(["%d" %(i+1) for i in s])
    #print "Exset * poo = %s;" % B
    return B

def get_random_set(end, incLen, exLen):
    site = 1
    poo = ""
    #poo = "exset * rand = "
    while site < end:
        s1 = site + random.randint(0,incLen)
        if s1 < end:
            s2 = s1 + random.randint(1,exLen)
            if s2 > end:
                s2 = end
            poo += "%s-%s " % (s1, s2)
            site = s2 + 1
        elif s1 == end:
            poo += "%s" % (s1)
            site = end
        else:
            site = end
    return poo

def get_random_set_triples(end, incLen, exLen):
    site = 1
    endC = end / 3
    poo = ""
    #poo = "exset * rand = "
    while site < endC:
        s1 = site + random.randint(0,incLen)
        if s1 < endC:
            s2 = s1 + random.randint(1,exLen)
            if s2 > endC:
                s2 = endC
            poo += "%s-%s " % (s1 * 3 - 2, s2 * 3)
            site = s2 + 1
        elif s1 == endC:
            poo += "%s-." % (s1)
            site = endC
        else:
            site = endC
    return poo      

def write_paupblock_for_random_data(numreps):
    for rep in range(0,numreps):
        print "incl all; restore all;"
        exclProp = random.uniform(EXCL_MIN,EXCL_MAX)
        print "exset * poo = " + get_better_random_set_triples(NCHAR, exclProp) + ";"
        inc = random.randint(1,10)
        ex = random.randint(1,100)
        exclProp = random.uniform(DEL_MIN,DEL_MAX)
        print "delete " + get_better_random_set(NTAX, exclProp) + ";"
        print "export file = rand.%03d.nex format=nex;" % (rep)
        
#start of actual code to execute

if type == "d":
    write_paupblock_for_random_data(num)
    sys.exit(0)

gen = GeneralSection()
gen.set_default()
#reset some to avoid too much insanity
gen.replace(TextOption("streefname", ["stepwise", "random"]))
gen.replace(TextOption("outputeachbettertopology", ["0"]))
gen.replace(TextOption("outputmostlyuselessfiles", ["0"]))

modN = ModelSection()
modN.set_default_nuc()

modC = ModelSection()
modC.set_default_codon()

modAA = ModelSection()
modAA.set_default_aa()

mast = MasterSection()
mast.set_default()
mast.replace(TextOption("intervallength", ["100"]))
mast.replace(IntOption("intervalstostore", 5, 5))

'''
print "[general]"
print gen.choose()
print modN.choose()
print "[master]"
print mast.choose()
'''

for i in range(0,num):
    filename = ("%03d" % i) + ".conf"
    f = open(filename, 'w')
    f.write("[general]\n")
    f.write("datafname = rand.%03d.nex\n" % i)
    gen.replace(TextOption("ofprefix", ["run.%03d" % i]))
    f.write(gen.choose())
    f.write("\n")
    if(type == "n"):
        f.write(modN.choose())
    elif(type == "c"):
        f.write(modC.choose())
    elif(type == "a"):
        f.write(modAA.choose())
    f.write("\n")
    f.write("[master]\n")
    f.write(mast.choose())
    f.close()
    
    '''
    f.write("[general]\n")
    #f.write(datafname.choose())
    f.write("datafname = rand.%03d.nex\n" % i)
    f.write(constraintfile.choose())
    f.write(streefname.choose())
    f.write(attachmentspertaxon.choose())
    #f.write(ofprefix.choose())
    f.write("ofprefix = run%03d\n" % i)
    f.write(randseed.choose())
    f.write(availablememory.choose())
    f.write(logevery.choose())
    f.write(saveevery.choose())
    f.write(refinestart.choose())
    f.write(outputeachbettertopology.choose())
    f.write(outputcurrentbesttopology.choose())
    f.write(enforcetermconditions.choose())
    f.write(genthreshfortopoterm.choose())
    f.write(scorethreshforterm.choose())
    f.write(significanttopochange.choose())
    f.write(outputphyliptree.choose())
    f.write(outputmostlyuselessfiles.choose())
    f.write(writecheckpoints.choose())
    f.write(restart.choose())
    f.write(outgroup.choose())
    f.write(searchreps.choose())
    f.write(collapsebranches.choose())
    f.write(outputsitelikelihoods.choose())
    f.write("\n")
    f.write(datatype.choose())
    f.write(ratematrix.choose())
    f.write(statefrequencies.choose())
    f.write(ratehetmodel.choose())
    f.write(numratecats.choose())
    f.write(invariantsites.choose())
    f.write("\n")
    f.write("[master]\n")
    f.write(nindivs.choose())
    f.write(holdover.choose())
    f.write(selectionintensity.choose())
    f.write(holdoverpenalty.choose())
    f.write(stopgen.choose())
    f.write(stoptime.choose())
    f.write("\n")
    f.write(startoptprec.choose())
    f.write(minoptprec.choose())
    f.write(numberofprecreductions.choose())
    f.write(treerejectionthreshold.choose())
    f.write(topoweight.choose())
    f.write(modweight.choose())
    f.write(brlenweight.choose())
    f.write(randnniweight.choose())
    f.write(randsprweight.choose())
    f.write(limsprweight.choose())
    f.write(intervallength.choose())
    f.write(intervalstostore.choose())
    f.write("\n")
    f.write(limsprrange.choose())
    f.write(meanbrlenmuts.choose())
    f.write(gammashapebrlen.choose())
    f.write(gammashapemodel.choose())
    f.write(uniqueswapbias.choose())
    f.write(distanceswapbias.choose())
    f.write("\n")
    f.write(bootstrapreps.choose())
    f.write(resampleproportion.choose())
    f.write(inferinternalstateprobs.choose())
    f.write("\n")
    f.close()
    '''

sys.exit(0)

print "[general]"
print datafname.choose()
print constraintfile.choose()
print streefname.choose()
print attachmentspertaxon.choose()
#print ofprefix.choose()
print "ofprefix = run%s" % i
print randseed.choose()
print availablememory.choose()
print logevery.choose()
print saveevery.choose()
print refinestart.choose()
print outputeachbettertopology.choose()
print outputcurrentbesttopology.choose()
print enforcetermconditions.choose()
print genthreshfortopoterm.choose()
print scorethreshforterm.choose()
print significanttopochange.choose()
print outputphyliptree.choose()
print outputmostlyuselessfiles.choose()
print writecheckpoints.choose()
print restart.choose()
print outgroup.choose()
print searchreps.choose()
print collapsebranches.choose()
print outputsitelikelihoods.choose()
print 
print datatype.choose()
print ratematrix.choose()
print statefrequencies.choose()
print ratehetmodel.choose()
print numratecats.choose()
print invariantsites.choose()
print
print "[master]"
print nindivs.choose()
print holdover.choose()
print selectionintensity.choose()
print holdoverpenalty.choose()
print stopgen.choose()
print stoptime.choose()
print
print startoptprec.choose()
print minoptprec.choose()
print numberofprecreductions.choose()
print treerejectionthreshold.choose()
print topoweight.choose()
print modweight.choose()
print brlenweight.choose()
print randnniweight.choose()
print randsprweight.choose()
print limsprweight.choose()
print intervallength.choose()
print intervalstostore.choose()
print
print limsprrange.choose()
print meanbrlenmuts.choose()
print gammashapebrlen.choose()
print gammashapemodel.choose()
print uniqueswapbias.choose()
print distanceswapbias.choose()
print
print bootstrapreps.choose()
print resampleproportion.choose()
print inferinternalstateprobs.choose()
f.close()
###

#OutputConfig()

quit()

if len(sys.argv) > 2:
    filename = sys.argv[1]
    file = open(filename, 'rU')
    perLine = sys.argv[2]
    if len(sys.argv) > 3:
        offset = int(sys.argv[3])
else:
    print "usage\nchopUpSeq.y <filename> <#per line> <optional codon NUMBER offset for first line>"
    sys.exit(0)



opt = ["streefname", ["stepwise", "random"]]
