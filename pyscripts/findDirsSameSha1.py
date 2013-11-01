#!/usr/bin/env python

import hashlib, sys, os
from os.path import normpath, walk, isdir, isfile, dirname, basename, \
    exists as path_exists, join as path_join

def recursive_superdir_generator(path):
    if not path:
        return
    yield path
    for sup in recursive_superdir_generator(os.path.split(path)[0]):
        yield sup


def is_subdir(path, superdir_set):
    if isinstance(superdir_set, list):
        superdir_set = set(superdir_set)
    if not superdir_set:
        return False
    shortest = min([len(superdir) for superdir in superdir_set])
    for sup in recursive_superdir_generator(path):
        if len(sup) < shortest:
            return False
        #print path, sup, superdir_set
        if sup in superdir_set:
            #print "subdir!", sup, path
            return True
    return False
        

def recursive_subdir_generator(path):
    contents = os.listdir(path)
    if not contents:
        return
    yield path
    for sub_dir in [ path_join(path, sub_obj) for sub_obj in contents if isdir(path_join(path, sub_obj)) and sub_obj[0] != '.' ]:
        for sub_sub_dir in recursive_subdir_generator(sub_dir):
            yield sub_sub_dir


def path_checksum(paths):
    """
    Recursively calculates a checksum representing the contents of all files
    found with a sequence of file and/or directory paths.

    """
    if not hasattr(paths, '__iter__'):
        raise TypeError('sequence or iterable expected not %r!' % type(paths))

    def _update_checksum(checksum, dirname, filenames):
        for filename in sorted(filenames):
            path = path_join(dirname, filename)
            if isfile(path):
                #print path
                fh = open(path, 'rb')
                while 1:
                    buf = fh.read(4096)
                    if not buf: 
                        break
                    checksum.update(buf)
                fh.close()

    chksum = hashlib.sha1()

    for path in sorted([normpath(f) for f in paths]):
        if path_exists(path):
            if isdir(path):
                walk(path, _update_checksum, chksum)
            elif isfile(path):
                _update_checksum(chksum, dirname(path), basename(path))

    return chksum.hexdigest()

hash_dict = {}
ignore_paths = set()
count = 0
skip_count = 0
match_count = 0
for top_dir in sys.argv[1:]:
    if isdir(top_dir):
        for dirname in recursive_subdir_generator(top_dir):
        #for dir in filtered_recursive_subdir_generator(top_dir, ignore_paths):
            #print ignore_paths, dirname 
            dirname = normpath(dirname)
            if not is_subdir(dirname, ignore_paths):
                checksum = path_checksum([dirname])
                #print checksum, dirname
                count += 1
                if checksum in hash_dict:
                    #don't consider a match between a directory and it's lone subdirectory
                    #print dirname, set(hash_dict[checksum])
                    if not is_subdir(dirname, hash_dict[checksum]):
                        hash_dict[checksum].append(dirname)
                        match_count += 1
                        #this will keep it from recursing any further if a higher level dir matched
                        ignore_paths.add(dirname)
                        print hash_dict[checksum]
                else:
                    hash_dict[checksum] = [dirname]
            else:
                skip_count += 1
            if count % 50 == 0:
                print count, skip_count, match_count

#print hash_dict
for sum, dirs in hash_dict.items():
    if len(dirs) > 1:
        print sum, dirs


