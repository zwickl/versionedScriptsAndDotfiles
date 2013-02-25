#!/usr/bin/env python
import sys
import re
import subprocess

#this just reads the specified files line by line, parses out the specified column, and 
#then pipes that to grep.  This allows full use of grep flags that are passed to this script.

arguments = []
patt = ''
files = []

try:
    column = int(sys.argv[1])
except ValueError:
    exit("first argument must be the column number to match the regex to")

foundPatt = False
for arg in sys.argv[2:]:
    if arg[0] == '-':
        arguments.append(arg)
    elif not patt:
        patt = arg
        foundPatt = True
    else:
        files.append(arg)

for f in files:
    ignored = 0
    for line in open(f, 'rb'):
        spline = line.strip().split()
        if len(spline) < column:
            ignored += 1
        else:
            el = re.escape(spline[column - 1])
            process = subprocess.Popen(' '.join(["echo", el, "|"] + ["grep"] + arguments + [ "'%s'" % patt ]), 
                    shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout_list, stderr_list = process.communicate()

            if stderr_list:
                print 'Child process exited with error? %s' % ''.join(stderr_list)
                exit()

            if stdout_list:
                sys.stdout.write(line)
    sys.stderr.write('ignored %d short lines\n' % ignored)

