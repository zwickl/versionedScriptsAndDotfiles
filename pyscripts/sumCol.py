#!/usr/bin/env python
import sys
from argparse import ArgumentParser

#use argparse module to parse commandline input
parser = ArgumentParser(description='sum or do other summaries of a column of numbers')

parser.add_argument("-s", "--sum", action="store_true", dest="outputSum", default=False, help="Output Sum")

parser.add_argument("--ignore-non-numeric", action="store_true", dest="ignoreNonNumeric", default=False, help="Ignore any column elements that can't be converted to floats")

parser.add_argument("-m", "--mean", action="store_true", dest="outputAve", default=False, help="Output Mean")

parser.add_argument("-r", "--range", action="store_true", dest="outputMinMax", default=False, help="Output Min and Max")

parser.add_argument("-a", "--all",  action="store_true", dest="outputAll", default=False, help="Output All Statistics")

parser.add_argument("-c", "--column", dest="columnNum", default=None, type=int, help="choose the column to output")

#variable number of arguments
parser.add_argument('filenames', nargs='*', default=[], 
                   help='a list of filenames to search (none for stdin)')

#now process the command line
options = parser.parse_args()

if options.outputAll:
    oSum = True
    oAve = True
    oMinMax = True
else:
    oSum = options.outputSum
    oAve = options.outputAve
    oMinMax = options.outputMinMax

    #if no options were entered, assume sum, otherwise -s is required to output the sum
    if not oAve and not oMinMax:
        oSum = True

if len(options.filenames) == 1:
    filename = options.filenames[0]
    try:
        infile = open(filename, 'rb')
    except IOError:
        raise RuntimeError("problem opening file")
elif options.filenames:
    sys.exit("expecting only one filename, got %s" % len(options.filenames))
else:
    infile = sys.stdin

lines = [ line.strip().split() for line in infile ]

if options.columnNum is None:
    if lines and len(lines[0]) > 1:
        raise RuntimeError('must pass either a single column or use the --column flag')
    else:
        options.columnNum = 1

colIndex = options.columnNum - 1
col = []
for line in lines:
    try:
        val = float(line[colIndex])
    except ValueError:
        if options.ignoreNonNumeric:
            sys.stderr.write('ignoring element \'%s\'\n' % line[colIndex])
            val = None
        else:
            raise
    except IndexError:
        exit('could not find column %d of line\n%s' % (options.colNum, line))
    #need to allow val to be zero here
    if val not in [ None, '' ]:
        col.append(val)

if not col:
    sys.stderr.write('No valid values read!\n')
    exit()

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
