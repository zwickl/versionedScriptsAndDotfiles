#!/usr/bin/env python
import os
import sys
import string
import StringIO

def usage():
	print "Usage: <script> filename [column to use] [histogram width]"

def makeHisto(lines, width = 1.0, col = 0):
	#figure out the range of the histogram
	if len(lines) is 0:
		sys.exit("list is empty")
	num = float(len(lines))
	if col > len(lines[0]):
		sys.exit("col # (%d) > # col in list (%d)" % (col, len(lines[0])))
	focalCol = [float(line[col]) for line in lines]
	colMin = min(focalCol)
	colMax   = max(focalCol)
	start = int(colMin - 1)
	while start + width < colMin:
		start = start + width
	end = int(colMax + 1)
	while end > colMax:
		end = end - width
	
	sys.stderr.write("range:\n%f %f %f %f\n" % (start, colMin, colMax, end+width))

#	counts = {}
	val = start
	bins = []
	bounds =[]
	while val <= (end + width * 0.5):
		#counts.append((val, 0))
		bins.append(0)
		bounds.append(val)
#		counts[val] = 0
		val = val + width
#		sys.stderr.write("%f " % val);
#	print counts
#	print bins
#	print bounds
	print "binStart\tbinEnd\tcount\tpercent"
	for element in focalCol:
		b = 0
	#	print element
		while bounds[b] + width < element:
			val = val + width
			b = b + 1
		bins[b] = bins[b] + 1
	for c in range(0, len(bins)):
		print "%f\t%f\t%d\t%f" % (bounds[c], bounds[c] + width, bins[c], float(bins[c]) / num)
'''
		val = start
		while val + width < element:
			val = val + width
		try:
			counts[val] = counts[val] + 1
		except KeyError:
			sys.exit("val %f not found in dict" % val )
			
	bins = counts.keys()
	bins.sort()

	for bin in bins:
		print "%f\t%d\t%f" % (bin, counts[bin], float(counts[bin]) / num)

'''
if len(sys.argv) < 2:
	usage()
	sys.exit(1)
fname = sys.argv[1]

width = 1.0
col = 0
try:
	col = int(sys.argv[2])
except:
	pass
try:
	width = float(sys.argv[3])
except:
	pass

try:
	file = open(fname, "ru")
except:
	sys.exit("problem opening file %s" % fname)

makeHisto([line.split() for line in file], width, col)
file.close()
