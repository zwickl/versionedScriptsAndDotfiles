#!/usr/bin/env python
import os
import sys
import string
import StringIO

from optparse import OptionParser

parser = OptionParser(add_help_option=True)

#parser.add_option("-h", "--help", dest="helpflag", help="Print this help" )
parser.add_option("-p", "--precision", dest="prec", type="int", default=2, help="Decimal precision of output")
parser.add_option("--horizontal", dest="hor", action="store_true", default=False, help="Separate cells with horizontal lines")
parser.add_option("--vertical", dest="vert", action="store_true",  default=False, help="Separate cells with vertical lines")
parser.add_option("--rowheaders", dest="rowfile", help="File containing row lables")
parser.add_option("--colheaders", dest="colfile", help="File containing column lables")

(options, args) = parser.parse_args()
#if options.helpflag != None:
#    print_usage("")

if len(args) == 1:
	filename = args[0]
	try:
		file = open(filename, 'r')
	except:
		sys.exit("problem opening file")
elif len(args) > 0:
	sys.exit("expecting only one argument, got %s" % args)
else:
	file = sys.stdin

outprec = int(options.prec)
haveRows = bool(options.rowfile is not None)
haveCols = bool(options.colfile is not None)
hor = options.hor
vert = options.vert
rowFile = None
colFile = None
rowHeads = []
colHeads = []

if haveRows:
	try:
		rowFile = open(options.rowfile, "rU")
	except:
		exit("Problem opening file %s" % options.rowfile)
	for r in rowFile:
		rowHeads.extend(string.replace(s, "_", " ") for s in r.split())
if haveCols:
	try:
		colFile = open(options.colfile, "rU")
	except:
		exit("Problem opening file %s" % options.colfile)
	for c in colFile:
		colHeads.extend(string.replace(s, "_", " ") for s in c.split())

lines = [ l.split() for l in file ]
numCols = len(lines[0])
floatRows = [[ float(e) for e in l ] for l in lines]
actualRows = [[ "%.*f" % (outprec, e) for e in l ] for l in floatRows ]

if haveRows == True:
	numCols = numCols + 1
	for l in range(0, len(actualRows)):
		actualRows[l].insert(0, rowHeads[l])
		#actualRows[l].insert(0, rowHeads[l].replace("_", " "))

print "\\begin {table}\n\\begin{center}\n\\resizebox{4in}{!}{"

if vert == True:
	print "\\begin{tabular}{ | " + " | ".join("c" for x in range(0, numCols)) + " | }"
else:
	print "\\begin{tabular}{ " + " ".join("c" for x in range(0, numCols)) + " }"

if hor == True:
	print "\\hline"
if haveCols == True:
	print " & ".join(colHeads) + "\\\\"
	print "\\hline"
elif hor == True:
	print "\\hline"

outRows = []
for l in actualRows:
	outRows.append(" & ".join(l) + "\\\\")
	if hor == True:
		outRows.append(" \\hline")
print "\n".join(outRows)
#for e in floatRows:
#	textLines.append(" & " + " & ".join( [ "%.*f" % (outprec, e) for e in l ] ))

#textLines = "\n".join([ " & " + " & ".join( [ "%.*f" % (outprec, e) for e in l ] ) + " \\\\" for l in floatRows ])
#print textLines
print "\\end{tabular}\n}\n\\end{center}\n\\end{table}"

