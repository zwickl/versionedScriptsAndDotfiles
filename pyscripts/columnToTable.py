#!/usr/bin/env python
import sys
import doctest

from argparse import ArgumentParser

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

def StringToIntOrFloat(string):
    """name says it all
    #>>> StringToIntOrFloat("1.2")
    #1.2
    """
    try:
        ret = int(string)
    except ValueError:
        ret = float(string)
    except:
        sys.stderr.write("string %s could not be interpreted as an int or float\n" % string)
        sys.exit(1)
    return ret

class row(object):
    
    def __init__(self, rowStr):
        #this saves the entries as strings
        self.colEntries = rowStr.strip().split()

    def element(self, c):
        if c >= len(self.colEntries):
            mess = 'requested index (%d) > than length of list (%d)\nline was %s\n' % (c, len(self.colEntries), str(self.colEntries))
            sys.stderr.write(mess)
            raise IndexError(mess)
        return self.colEntries[c]

    def columns(self):
        return len(self.colEntries)

    def __print__(self):
        print '\t'.join(self.colEntries)
        #for ce in self.colEntries:
        #    print ce,

class rowset(object):
    #rows = []
    def __init__(self):
        self.rows = []

    def append_row(self, rStr):
        self.rows.append(row(rStr))

    def column_to_list(self, cnum, missingString=None):
        clist = []
        for poo in self.rows:
            try:
                clist.append(poo.element(cnum))
            except IndexError:
                if missingString:
                    clist.append(missingString)
                else:
                    raise
        #clist = [ poo.element(cnum) for poo in self.rows ]
        return clist

    def columns(self):
        return self.rows[0].columns()

    def clear(self):
        self.rows = []

    def __print__(self):
        for r in self.rows:
            r.__print__()
            #print

#this just runs the tests embedded in doc strings
#if __name__ == "__main__":
#    doctest.testmod()

parser = ArgumentParser()

parser.add_argument("-i", "--input-file", dest="filename", default=None, type=str, help="input file")
parser.add_argument("-m", "--missing-string", dest="missingString", default=None, type=str, 
        help="allow missing values, and replace with specified string")
parser.add_argument("-r", "--rows", dest="nrows", type=int, required=True, help="number of rows per chunk")
parser.add_argument("-f", "--first-cols", dest="firstCols", nargs="*", type=int, default=None,
        help="columns to use in first row set")
parser.add_argument("-c", "--cols", dest="cols", nargs="*", type=int, default=None, help="columns to use in successive row sets")
parser.add_argument("-t", "--transpose", action="store_true", dest="transflag", default=False, help="transpose matrix")

options = parser.parse_args()

if options.filename is None:
    sys.stderr.write("No filename passed (-i), assuming stdin\n")

infile = open(options.filename, "ri") if options.filename else sys.stdin

chunks = []
rset = rowset()
rnum = 0
for line in infile:
    rset.append_row(line)
    rnum = rnum + 1
    if rnum == options.nrows:
        chunks.append(rset)
        rset = rowset()
        rnum = 0

#if specific columns weren't passed in, use them all.
#remember that it assumes that they are specified starting at 1
if options.firstCols is None:
    firstCols = range(1, chunks[0].columns() + 1)
else:
    firstCols = options.firstCols

if options.cols is None:
    cols = range(1, chunks[1].columns() + 1)
else:
    cols = options.cols

sys.stderr.write("Found %d total sets\n" % len(chunks))
sys.stderr.write("using columns %s from the first set\n" % ", ".join([str(f) for f in firstCols]))
sys.stderr.write("and columns %s from later sets\n" % ", ".join([str(f) for f in cols]))

#get the proper columns from the first set
tableCols = [chunks[0].column_to_list(c - 1, options.missingString) for c in firstCols]

for thisChunk in chunks[1:]:
    tableCols.extend(thisChunk.column_to_list(c - 1, options.missingString) for c in cols)

if not options.transflag:
    #nested list comp
    print "\n".join(["\t".join([col[r] for col in tableCols]) for r in xrange(options.nrows)])
    
else:
    #nested list comp
    print "\n".join(["\t".join([col[r] for r in xrange(options.nrows)]) for col in tableCols])

