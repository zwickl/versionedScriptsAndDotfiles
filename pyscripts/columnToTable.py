#!/usr/bin/env python
import sys

from argparse import ArgumentParser, FileType

def print_usage(errMess):
    if errMess:
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
    sys.exit()


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


class rowset(object):
    def __init__(self):
        self.rows = []

    def append_row(self, rStr):
        self.rows.append(row(rStr))

    def column_to_list(self, cnum, missing_string=None):
        clist = []
        for poo in self.rows:
            try:
                clist.append(poo.element(cnum))
            except IndexError:
                if missing_string:
                    clist.append(missing_string)
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

    def __len__(self):
        return len(self.rows)

parser = ArgumentParser()

parser.add_argument("-i", "--input", default=sys.stdin, type=FileType('rU'), help="input file")
parser.add_argument("-m", "--missing-string", default=None, type=str, 
                    help="allow missing values, and replace with specified string")

parser.add_argument("-r", "--rows", dest="nrows", type=int, required=True, help="number of rows per chunk")

parser.add_argument("-f", "--first-cols", nargs="*", type=int, default=None,
                    help="columns to use in first row set")

parser.add_argument("-c", "--cols", nargs="*", type=int, default=None, help="columns to use in successive row sets")

parser.add_argument("-t", "--transpose", action="store_true", default=False, help="transpose matrix")
parser.add_argument("-q", "--quiet", action="store_true", dest="quiet", default=False, help="no messages output to the screen")
parser.add_argument("--skip-rows", type=int, default=0, help="initial rows to skip (NOTE: the --rows specification should not include these rows)")

options = parser.parse_args()

chunks = []
rset = rowset()
for num, line in enumerate(options.input):
    if num >= options.skip_rows:
        rset.append_row(line)
        if len(rset) == options.nrows:
            chunks.append(rset)
            rset = rowset()

if rset:
    if not options.quiet:
        sys.stderr.write('WARNING - last set of rows has fewer than --rows value.  Check --rows\n')
    chunks.append(rset)

#if specific columns weren't passed in, use them all.
#remember that it assumes that they are specified starting at 1
firstCols = options.first_cols or range(1, chunks[0].columns() + 1)
cols = options.cols or range(1, chunks[1].columns() + 1)

if not options.quiet:
    sys.stderr.write("Found %d total sets\n" % len(chunks))
    sys.stderr.write("using columns %s from the first set\n" % ", ".join([str(f) for f in firstCols]))
    sys.stderr.write("and columns %s from later sets\n" % ", ".join([str(f) for f in cols]))

#get the proper columns from the first set
tableCols = [chunks[0].column_to_list(c - 1, options.missing_string) for c in firstCols]

if len(chunks) > 1:
    for thisChunk in chunks[1:]:
        tableCols.extend(thisChunk.column_to_list(c - 1, options.missing_string) for c in cols)

if not options.transpose:
    #nested list comp
    print "\n".join(["\t".join([col[r] for col in tableCols]) for r in xrange(options.nrows)])
    
else:
    #nested list comp
    print "\n".join(["\t".join([col[r] for r in xrange(options.nrows)]) for col in tableCols])

