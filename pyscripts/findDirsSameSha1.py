#!/usr/bin/env python

import hashlib, sys, os
from os.path import normpath, walk, isdir, isfile, dirname, basename, \
    exists as path_exists, join as path_join
from itertools import chain

def recursive_superdir_generator(path):
    if not path:
        return
    yield path
    for super in recursive_superdir_generator(os.path.split(path)[0]):
        yield super


def is_subdir(path, superdir_set):
    #if isinstance(superdir_list, str):
    #    superdir_list = [superdir_list]
    for super in recursive_superdir_generator(path):
        if super in superdir_set:
            print "subdir!", super, path
            return True
    return False
        
    exit()
    
    for super in superdir_list:
        if super + '/' in path:
            return True
    return False


def recursive_subdir_generator(path):
    yield path
    for sub_dir in [ path_join(path, sub_obj) for sub_obj in os.listdir(path) if isdir(path_join(path, sub_obj)) and sub_obj[0] != '.' ]:
        for sub_sub_dir in recursive_subdir_generator(sub_dir):
            yield sub_sub_dir


def filtered_recursive_subdir_generator(path, skips):
    print path, skips
    for skip in skips:
        if skip in path:
            return
    yield path
    for sub_dir in [ path_join(path, sub_obj) for sub_obj in os.listdir(path) if isdir(path_join(path, sub_obj)) and sub_obj[0] != '.' ]:
        for sub_sub_dir in [ ssd for ssd in filtered_recursive_subdir_generator(sub_dir, skips) if ssd is not None ]:
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
        for dir in recursive_subdir_generator(top_dir):
        #for dir in filtered_recursive_subdir_generator(top_dir, ignore_paths):
            if not is_subdir(dir, ignore_paths):
                checksum = path_checksum([dir])
                count += 1
                if checksum in hash_dict:
                    #don't consider a match between a directory and it's lone subdirectory
                    ok = True
                    for found_dir in hash_dict[checksum]:
                        if found_dir + '/' in dir:
                            ok = False
                            break
                    if ok:
                        hash_dict[checksum].append(dir)
                        match_count += 1
                        #this will keep it from recursing any further if a higher level dir matched
                        ignore_paths.add(dir)
                        #print hash_dict[checksum]
                else:
                    hash_dict[checksum] = [dir]
            else:
                skip_count += 1
            if count % 10 == 0:
                print count, skip_count, match_count

#print hash_dict
for sum, dirs in hash_dict.items():
    if len(dirs) > 1:
        print sum, dirs


