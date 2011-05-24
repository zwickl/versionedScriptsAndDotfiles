#!/usr/bin/env python
import os
import string

def islogfile(x): return x.endswith('.log00.log')

def numeric_compare(x,y):
	f=x[0:2]
	s=y[0:2]
	return int(f)-int(s)

runnum = 1
files = os.listdir('.')
files = filter(islogfile, files)
files.sort(cmp=numeric_compare)
if len(files) < 1:
	print "Sorry, no log files found"

file = open('expected.scr', 'r')
expected = file.readlines()
file.close()

num=0
for i in files:
	if num >= len(expected): break
	if i[0].isdigit():
		scores = []
		file = open(i, 'r')
		line = file.readline()
		while line != "":
			if line.startswith('Final'):
				divided = line.split()
				scores.append(string.atof(divided[1]))
			line = file.readline()
#		print scores
	#	sum = string.atof(expected[num]) + string.atof(divided[1])
		if(len(scores) > 0):
			sum = string.atof(expected[num]) + max(scores)
			#print '%(f)s %(e)f %(g)s' % {'f': i, "e": string.atof(expected[num]), "g": divided[1]}
			if(abs(sum) < 0.05):
				print '%(f)25s %(s)15f PASS ' % {'f': i, "s": sum}, 
			else:
				print '%(f)25s %(s)15f FAIL' % {'f': i, "s": sum},
			print expected[num] , scores
		else:
			print i, "no score found in file?"
		file.close()
		num = num + 1






