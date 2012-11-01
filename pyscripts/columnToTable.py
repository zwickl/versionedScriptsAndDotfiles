#!/usr/bin/env python
import os
import sys
import string
import StringIO
import doctest

from optparse import OptionParser

#args are the name of the option, the actual flag that was used, the value, and the parser itself
#Note that I'm not sure if assigning the split back to value is a good practice, but it works
def comma_split(option, opt_str, value, parser):
    intList = []
    strList = value.split(",")
    for s in strList:
        intList.append(int(s))
        setattr(parser.values, option.dest, intList)

def print_usage(errMess):
    if len(errMess) != 0:
        sys.stderr.write(errMess)
    parser.print_help()
    sys.stderr.write("i.e., given this inp\n")
    sys.stderr.write("1    a    poo\n")
    sys.stderr.write("2    b    moo\n")
    sys.stderr.write("3    c    doo\n")
    sys.stderr.write("1    a    roo\n")
    sys.stderr.write("2    b    goo\n")
    sys.stderr.write("3    c    noo\n")
    sys.stderr.write("this: columnToTable.py -i filename -r 3 -f 1,2,3 -c 3\n")
    sys.stderr.write("would give:\n")
    sys.stderr.write("1    a    poo    roo\n")
    sys.stderr.write("2    b    moo    goo\n")
    sys.stderr.write("3    c    doo    noo\n")
    exit()

def StringToIntOrFloat(str):
    """name says it all
    #>>> StringToIntOrFloat("1.2")
    #1.2
    """
    try:
        ret = int(str)
    except ValueError:
        ret = float(str)
    except:
        sys.stderr.write("string %s could not be interpreted as an int or float\n" % str)
        sys.exit(1)
    return ret

class row:
    
    def __init__(self, rowStr):
        #this saves the entries as strings
        self.colEntries = [x for x in rowStr.split()]
        
        #this would, obviously, always force floats
        #self.colEntries = [float(x) for x in rowStr.split()]
        
        #this converts to ints or floats at necessary
        #self.colEntries = [StringToIntOrFloat(x) for x in rowStr.split()]
        #print [type(x) for x in self.colEntries]

    def element(self, c):
        if c >= len(self.colEntries):
            mess = 'requested index (%d) > than length of list (%d)' % (c, len(self.colEntries))
            sys.stderr.write(mess)
            raise IndexError(mess)
        return self.colEntries[c]

    def columns(self):
        return len(self.colEntries)

    def __print__(self):
        for ce in self.colEntries:
            print ce,

class rowset:
    #rows = []
    def __init__(self):
        self.rows = []

    def append_row(self, rStr):
        self.rows.append(row(rStr))

    def column_to_list(self, cnum):
        clist = [ poo.element(cnum) for poo in self.rows ]
        return clist

    def columns(self):
        return self.rows[0].columns()

    def clear(self):
        self.rows = []

    def __print__(self):
        for r in self.rows:
            r.__print__()
            print

parser = OptionParser(add_help_option=False)

parser.add_option("-h", "--help", dest="helpflag")
parser.add_option("-i", "--inputfile", dest="filename", type="string", help="input file", metavar="FILE")
parser.add_option("-r", "--rows", dest="nrows", type="int", help="number of rows per chunk", metavar="#")
parser.add_option("-f", "--firstcols", dest="firstCols", type="string", action="callback", callback=comma_split, help="columns to use in first row set", metavar="# # ...")
parser.add_option("-c", "--cols", dest="cols", type="string", action="callback", callback=comma_split, help="columns to use in successive row sets", metavar="# # ...")
parser.add_option("-t", "--transpose", action="store_true", dest="transflag", default=False, help="transpose matrix")

(options, args) = parser.parse_args()
if options.helpflag != None:
    print_usage("")

#this just runs the tests embedded in doc strings
#if __name__ == "__main__":
#    doctest.testmod()

if options.filename == None:
    #print_usage("you must pass the -i option")
    sys.stderr.write("No filename passed (-i), assuming stdin\n")
if options.nrows == None:
    print_usage("you must pass the -r option")
'''
if options.firstCols == None:
    print_usage("you must pass the -f option")
if options.cols == None:
    print_usage("you must pass the -c option")
'''

transpose = options.transflag

infile = sys.stdin if options.filename == None else open(options.filename, "ri")

nrows = options.nrows
try:
    firstCols = options.firstCols
except:
    firstCols = -1

try:
    cols = options.cols
except:
    cols = -1

chunks = []
rset = rowset()
rnum = 0
for line in infile.readlines():
    rset.append_row(line)
    rnum = rnum + 1
    if rnum == nrows:
        chunks.append(rset)
        rset = rowset()
        rnum = 0


#if specific columns weren't passed in, use them all.
#remember that it assumes that they are specified starting at 1
if firstCols == None:
    firstCols = range(1, chunks[0].columns() + 1)
if cols == None:
    cols = range(1, chunks[1].columns() + 1)

sys.stderr.write("Found %d total sets\n" % len(chunks))
sys.stderr.write("using columns %s from the first set\n" % ", ".join([str(f) for f in firstCols]))
sys.stderr.write("and columns %s from later sets\n" % ", ".join([str(f) for f in cols]))

#get the proper columns from the first set
tableCols = [chunks[0].column_to_list(c - 1) for c in firstCols]
#print tableCols
#get the rest of the columns.  this is using a list comprehension for the side effects, not the list return value
#[ [tableCols.append(thisChunk.column_to_list(c - 1)) for c in cols] for thisChunk in chunks[1:] ]

for thisChunk in chunks[1:]:
    tableCols.extend([thisChunk.column_to_list(c - 1) for c in cols])
#print tableCols

if not transpose:
    #nested list comp
    print "\n".join(["\t".join([col[r] for col in tableCols]) for r in range(0,nrows)])
    
    '''one list comp
    for r in range(0,nrows):
    print "\t".join([col[r] for col in tableCols])
    '''
    
    '''plain loops
    for r in range(0,nrows):
    for col in tableCols:
        print col[r] + "\t",
    print
    '''
else:
    #nested list comp
    print "\n".join(["\t".join([col[r] for r in range(0,nrows)]) for col in tableCols])

    '''one list comp
    for col in tableCols:
    print "\t".join([col[r] for r in range(0,nrows)])
    '''

    '''plain loops
    for col in tableCols:
    for r in range(0,nrows):
        print "%f\t" % col[r],
    print
    '''
