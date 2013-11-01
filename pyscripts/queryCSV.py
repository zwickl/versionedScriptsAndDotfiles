#!/usr/bin/env python
import sys
import re
from csv import DictReader
from argparse import ArgumentParser

parser = ArgumentParser(description='extract information from a csv (or tsv) file based on a query')

parser.add_argument('-l', '--list', action='store_true', default=False,
                    help='list the fields (column headers) available for query or output')

parser.add_argument('-q', '--query', type=str, default='lambda row: row',
        help='Lambda to select rows. Should return True/False!, e.g. \'lambda row: row["treatment"].lower() == "mafft.cds"\' (default is to include all rows)')

parser.add_argument('-o', '--output', type=str, default=None,
        help='lambda to select which columns of each row to output, i.e. \'lambda row: [row["treatment"], row["majCount"], row["majSupp"]]\'')

parser.add_argument('filename', type=str, 
                    help='file to read from, required to have column headers')

opt = parser.parse_args()

lines = [line for line in open(opt.filename, 'rU')]
header_row = lines[0]

#remove any extra times the header row repeats in the file
lines = [header_row] + [ line for line in lines if line != header_row ]
header_row = header_row.split()

fileData = DictReader(lines, delimiter='\t')

if opt.list:
    print fileData.fieldnames
    sys.exit()

#filter csv row dicts based on some criterion
qfunc = eval(opt.query)
filtered = filter(qfunc, fileData)

if not opt.output:
    #if no specific output selection function was passed in, output everything in same column order
    print '\t'.join(header_row)
    for row in filtered:
        sys.stdout.write('%s\n' % '\t'.join([row[field] for field in header_row]))
else:
    ofunc = eval(opt.output)
    to_output = [ ofunc(line) for line in filtered ]
    for row in to_output:
        if isinstance(row, str):
            sys.stdout.write('%s\n' % row)
        else:
            sys.stdout.write('%s\n' % '\t'.join(row))

exit()


func = lambda row:  row['triplet'] == 'barthii.punctata.officinalis.brachyantha.dat' 

if len(sys.argv) > 2:
    func = lambda row: re.match(sys.argv[2], row['triplet']) is not None
else:
    func = lambda row: re.match('sat.*punc.*rufi.*sat', row['triplet']) is not None

#func = lambda row: [ row for row in fileData if row['triplet'] == 'barthii.punctata.officinalis.brachyantha.dat' ]
out = [ func(row) for row in fileData ]
out = filter(func, fileData)
#out = [ (row['treatment'], row['nameMaj']) for row in fileData if func(row) ]
out = [ row['nameMaj'] for row in fileData if func(row) ]

for o in out:
    print o

triplets = set([ row['triplet'] for row in fileData ])
for t in triplets:
    out = set([ row['nameMaj'] for row in fileData if row['triplet'] == t ])
    if len(out) > 1:
        print t
