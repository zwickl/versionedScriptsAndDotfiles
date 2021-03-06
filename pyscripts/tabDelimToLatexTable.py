#!/usr/bin/env python
import sys
import string

from optparse import OptionParser

parser = OptionParser(add_help_option=True)

#parser.add_option("-h", "--help", dest="helpflag", help="Print this help" )
parser.add_option("-p", "--precision", dest="prec", type="int", default=2, help="Decimal precision of output")
parser.add_option("--horizontal", dest="hor", action="store_true", default=False, help="Separate cells with horizontal lines")
parser.add_option("--vertical", dest="vert", action="store_true",  default=False, help="Separate cells with vertical lines")
parser.add_option("--rowheaders", dest="rowfile", default=None, help="File containing row labels")
parser.add_option("--colheaders", dest="colfile", default=None, help="File containing column labels")

(options, args) = parser.parse_args()

if len(args) == 1:
    filename = args[0]
    try:
        infile = open(filename, 'r')
    except IOError:
        sys.exit("problem opening file")
elif args:
    sys.exit("expecting only one argument, got %s" % args)
else:
    infile = sys.stdin

outprec = int(options.prec)
rowHeads, colHeads = [], []

if options.rowfile:
    try:
        rowFile = open(options.rowfile, "rU")
    except IOError:
        exit("Problem opening file %s" % options.rowfile)
    for r in rowFile:
        rowHeads.extend(string.replace(s, "_", " ") for s in r.split())
if options.colfile:
    try:
        colFile = open(options.colfile, "rU")
    except IOError:
        exit("Problem opening file %s" % options.colfile)
    for c in colFile:
        colHeads.extend(string.replace(s, "_", " ") for s in c.split())

lines = [ l.split() for l in infile ]
numCols = len(lines[0])
floatRows = [[ float(e) for e in l ] for l in lines]
actualRows = [[ "%.*f" % (outprec, e) for e in l ] for l in floatRows ]

if options.rowfile:
    numCols = numCols + 1
    for l in range(0, len(actualRows)):
        actualRows[l].insert(0, rowHeads[l])
        #actualRows[l].insert(0, rowHeads[l].replace("_", " "))

print "\\begin {table}\n\\begin{center}\n\\resizebox{4in}{!}{"

if options.vert:
    print "\\begin{tabular}{ | " + " | ".join("c" for x in range(0, numCols)) + " | }"
else:
    print "\\begin{tabular}{ " + " ".join("c" for x in range(0, numCols)) + " }"

if options.hor:
    print "\\hline"
if options.colfile:
    print " & ".join(colHeads) + "\\\\"
    print "\\hline"
elif options.hor:
    print "\\hline"

outRows = []
for l in actualRows:
    outRows.append(" & ".join(l) + "\\\\")
    if options.hor:
        outRows.append(" \\hline")
print "\n".join(outRows)
#for e in floatRows:
#    textLines.append(" & " + " & ".join( [ "%.*f" % (outprec, e) for e in l ] ))

#textLines = "\n".join([ " & " + " & ".join( [ "%.*f" % (outprec, e) for e in l ] ) + " \\\\" for l in floatRows ])
#print textLines
print "\\end{tabular}\n}\n\\end{center}\n\\end{table}"

