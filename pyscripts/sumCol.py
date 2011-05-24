#!/usr/bin/env python
import os
import sys
import string
import StringIO

from optparse import OptionParser

parser = OptionParser(add_help_option=False)

parser.add_option("-h", "--help", dest="helpflag")
parser.add_option("-s", "--sum", action="store_true", dest="outputSum", default=False, help="Output Sum")
parser.add_option("-m", "--mean", action="store_true", dest="outputAve", default=False, help="Output Mean")
parser.add_option("-r", "--range", action="store_true", dest="outputMinMax", default=False, help="Output Min and Max")
parser.add_option("-a", "--all",  action="store_true", dest="outputAll", default=False, help="Output All Statistics")

(options, args) = parser.parse_args()
if options.helpflag != None:
    print_usage("")

if options.outputAll is True:
	oSum = True
	oAve = True
	oMinMax = True
else:
	oSum = options.outputSum
	oAve = options.outputAve
	oMinMax = options.outputMinMax

	#if no options were entered, assume sum, otherwise -s is required to output the sum
	if oAve is False and oMinMax is False:
		oSum = True

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
'''
if len(sys.argv) > 1:
	filename = sys.argv[1]
	file = open(filename, 'r')
else:
	file = sys.stdin
'''
col = [float(line) for line in file]

cSum = sum(col)
cNum = len(col)

out = []

if oSum:
	out.append(str(cSum))
if oAve:
	out.append(str(cSum / cNum))
if oMinMax:
	out.append(str(min(col)))
	out.append(str(max(col)))

print "\t".join(out)
