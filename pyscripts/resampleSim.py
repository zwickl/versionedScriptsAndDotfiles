#!/usr/bin/env python
import random
n=500
data = range(n)
nreps = 1000
nhits = 0
for reps in range(nreps):
    data = range(n)
    nrounds = 1
    while True:
	homo = True
	x = data[0]
	for i in data:
		if i != x:
			homo = False
			break
	if homo:
		break	
	nrounds += 1
	data = [random.choice(data) for j in range(n)]
    if data[0] < n/2:
	nhits += 1
    print("%d out of %d in %d iters" % ( nhits, reps, nrounds))
