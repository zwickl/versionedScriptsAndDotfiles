#!/usr/bin/env python

for inc in range(1, 16):
    nums = range(1, 16)
    nums.remove(inc)
    print "exset * poo = " ,
    for n in nums:
        print n ,
    print ";"
    #[print n for n in nums]
 
