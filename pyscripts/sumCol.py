#!/usr/bin/env python
import sys
from argparse import ArgumentParser, FileType, RawDescriptionHelpFormatter

#use argparse module to parse commandline input
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter, description='''sum or do other summaries of columns of numbers\nexample: to sum PIDs of user processes listed by shell ps command:\nps | sumCol.py --column 1 --ignore-non-numeric''')

parser.add_argument("-s", "--sum", action="store_true", default=False, help="Output Sum")

parser.add_argument("-i", "--ignore-non-numeric", action="store_true", default=False, help="Ignore any column elements that can't be converted to floats")

parser.add_argument("-m", "--mean", action="store_true", default=False, help="Output Mean")

parser.add_argument("-r", "--range", action="store_true", default=False, help="Output Min and Max")

parser.add_argument("-a", "--all",  action="store_true", default=False, help="Output All Statistics")

parser.add_argument("-c", "--column", default=None, type=int, help="choose the column number to output (starting at 1)")

parser.add_argument("-p", "--precision", default=4, type=int, help="number of digits past decimal for floating point (default 4)")

parser.add_argument('infile', nargs='?', default=sys.stdin, type=FileType('rU'), help='filename to search (none for stdin)')

#now process the command line
options = parser.parse_args()

if options.all:
    oSum, oAve, oMinMax = True, True, True
else:
    oAve = options.mean
    oMinMax = options.range
    #if no options were entered, assume sum, otherwise -s is required to output the sum
    oSum = True if not (oAve or oMinMax) else options.sum

lines = [ line.strip().split() for line in options.infile ]

if options.column is None:
    if lines and len(lines[0]) > 1:
        raise RuntimeError('must pass either a single column or use the --column flag')
    else:
        options.column = 1

colIndex = options.column - 1
col = []
for line in lines:
    try:
        val = float(line[colIndex])
    except ValueError:
        if options.ignore_non_numeric:
            sys.stderr.write('ignoring element \'%s\'\n' % line[colIndex])
            val = None
        else:
            raise
    except IndexError:
        sys.exit('could not find column %d of line\n%s' % (options.column, line))
    #need to allow val to be zero here
    if val not in [ None, '' ]:
        col.append(val)

'''no lines read shouldn't really be a problem
if not col:
    sys.stderr.write('No valid values read!\n')
    exit()
'''

cSum = sum(col)
cNum = len(col)

if not cNum:
    if options.infile is sys.stdin:
        sys.exit('No values read from stdin!')
    else:
        sys.exit('No values read from file %s!' % options.infile)

else:
    out = []

    if oSum:
        out.append(str(cSum))
    if oAve:
        precString = ('%%.%df' % options.precision)
        out.append(precString % (cSum / cNum))
    if oMinMax:
        out.append(str(min(col)))
        out.append(str(max(col)))

    print(("\t".join(out)))
    
