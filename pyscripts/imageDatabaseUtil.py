#!/usr/bin/env python
import sys
import os
import time
import re
from datetime import date
import shelve
from hashlib import md5
import exifread
import argparse
import sqlite3
from dzutils import flatten_array


def md5Checksum(filePath):
    with open(filePath, 'rb') as fh:
        m = md5()
        while True:
            data = fh.read(8192)
            if not data:
                break
            m.update(data)
        return m.hexdigest()


class FileInfo(object):
    def __init__(self, fname):
        self.path = fname
        self.fname = os.path.basename(fname)
        self.dir = fname.split('/')[-2]
        #self.stat = os.stat(fname)
        #self.file_time = time.ctime(os.path.getctime(fname))
        self.file_time = str(date.fromtimestamp(os.path.getmtime(fname)))
        self.md5 = md5Checksum(fname)
        self.hash = hash(self.md5)

    def __hash__(self):
        #return self.md5
        return self.hash

    def __eq__(self, other):
        return self.md5 == other.md5

    def __lt__(self, other):
        return self.fname < other.fname


class ImageFile(FileInfo):
    def __init__(self, fname):
        super(ImageFile, self).__init__(fname)
        with open(fname, 'rb') as f:
            self.tags = exifread.process_file(f)
        if not 'EXIF DateTimeOriginal' in self.tags:
            #sys.stderr.write('No exif data for %s\n' % fname)
            self.exif_time = None
        else:
            match = re.search('[0-9]+:[0-9]+:[0-9]+', '%s' % self.tags['EXIF DateTimeOriginal'])
            if match:
                self.exif_time = re.sub(':', '-', match.group(0).strip())
            else:
                sys.stderr.write('Could not parse exif data %s: %s\n' % (fname, self.tags['EXIF DateTimeOriginal']))
                self.exif_time = None

    def __repr__(self):
        ret = '%50s\t%12s\t%12s\t%12s' % (self.path, self.exif_time, self.file_time, self.dir)
        if self.exif_time is None:
            ret += ' missing_exif '
        else:
            if self.exif_time != self.file_time:
                ret += ' exif_file_mis '
            if self.exif_time != self.dir:
                ret += ' exif_dir_mis '
        return ret

    def tuple_repr(self):
        ex_time = self.exif_time or 'None'
        #print self.exif_time
        return (self.path, ex_time, self.file_time, self.dir, self.md5)

def open_image_shelf(dbfile):
    if os.path.exists(dbfile):
        sys.stderr.write('opening existing database file %s\n' % dbfile)
        #img_dict = shelve.open(dbfile, protocol=-1, writeback=use_writeback)
        #known_paths = shelve.open(pfile, protocol=-1, writeback=use_writeback)
        
        img_shelf = shelve.open(dbfile)
        if 'known_paths' in img_shelf:
            known_paths = img_shelf['known_paths']
        if 'img_dict' in img_shelf:
            img_dict = img_shelf['img_dict']

        #sys.stderr.write('%d paths and %d unique files loaded\n' % (len(known_paths), len(img_dict)))

    else:
        sys.stderr.write('creating new database file %s\n' % dbfile)
        #img_dict = shelve.open(dbfile, protocol=-1, writeback=use_writeback)
        #known_paths = shelve.open(pfile, protocol=-1, writeback=use_writeback)
        
        img_shelf = shelve.open(dbfile, writeback=True)
        known_paths = set()
        img_dict = {}
        img_shelf['img_dict'] = img_dict
        img_shelf['known_paths'] = known_paths

    return (img_shelf, img_dict, known_paths)



parser = argparse.ArgumentParser(description='Create a database and manipulate image files')

parser.add_argument('-d', '--duplicates', action='store_true', default=False, help='output list of duplicate files to stdout')

parser.add_argument('-v', '--reverse', action='store_true', default=False, help='reverse sense of no-exif or mismatch match (files NOT meeting condition)')

parser.add_argument('--no-exif', action='store_true', default=False, help='output list of files without exif information to stdout')

parser.add_argument('--exif-file-mismatch', action='store_true', default=False, help='output list of files with conflict between file date and exif info to stdout')

parser.add_argument('--exif-dir-mismatch', action='store_true', default=False, help='output list of files with conflict between directory name and exif info to stdout')

parser.add_argument('--read-only', action='store_true', default=False, help='perform actions on the database (even adding files), but don\'t write it back to file (faster)')

parser.add_argument('--minimal', action='store_true', default=False, help='only output a list of filenames for the chosen option, no other information')

#parser.add_argument('-db', '--database', default='database.sql', help='name of database file to load/create (default database.shelf)')
parser.add_argument('-db', '--database', default='images', help='name of database file to load/create (default database.shelf)')

#parser.add_argument('-db2', '--database2', default='database2.sql', help='name of 2nd database file to load for set comparisons (default database2.shelf)')
parser.add_argument('-db2', '--database2', default='images2', help='name of 2nd database file to load for set comparisons (default database2.shelf)')

parser.add_argument('--action', default=None, choices=['intersection', 'union', 'xor', 'minus', 'rminus'], help='set operations to perform on two databases')

parser.add_argument('-l', '--list', default=None, help='list entries containing this string in the path field')

parser.add_argument('files', nargs='*', default=[], help='image or movie files to process')

options = parser.parse_args()

use_writeback = False

if not options.action:
    #img_shelf, img_dict, known_paths = open_image_shelf(options.database)

    conn = sqlite3.connect('database.sql')
    conn.row_factory = sqlite3.Row
    c = conn.cursor()
    c.execute('''CREATE TABLE if not exists %s
            (path TEXT PRIMARY KEY, exif_time TEXT, file_time TEXT, dir TEXT, md5 TEXT)''' % options.database)

    img_dict = {}
    
    numAdded = 0
    numHad = 0
    for fn in options.files:
        if os.path.isfile(fn):
                
            c.execute('SELECT * FROM %s WHERE path=?' % options.database, (fn,))
            if not c.fetchone():
                if numAdded == 0:
                    sys.stderr.write('adding files ...\n')
                numAdded += 1
                f = ImageFile(fn)
                img_dict.setdefault(f.md5, []).append(f)
                #known_paths.add(fn)
                if numAdded % 100 == 0:
                    sys.stderr.write('%d ' % numAdded)
            else:
                numHad += 1
        else:
            sys.stderr.write('%s is not a file. Skipping.' % fn)

    sys.stderr.write('\n%d items added to database\n%d items were already in database\n%d total items in database\n' % (numAdded, numHad, len(img_dict)))
    c.executemany("INSERT OR IGNORE INTO %s VALUES (?, ?, ?, ?, ?)" % options.database, [ img.tuple_repr() for img in flatten_array(img_dict.values()) ])
 
    if options.duplicates:
        c.execute("SELECT md5 FROM %s GROUP BY md5 HAVING count(md5) > 1" % options.database)
        for md5 in c.fetchall():
            #for row in c.execute("SELECT * FROM %s WHERE md5=?" % options.database, (md5[0],)):
            #    print row
            #print
            dupes = [ row for row in c.execute("SELECT * FROM %s WHERE md5=?" % options.database, (md5[0],)) ]
            #print '\n'.join([str(d) for d in dupes])
            for n1, d1 in enumerate(dupes):
                for d2 in dupes[n1 + 1:]:
                    if list(d1)[1:] == list(d2)[1:]:
                        print 'exact dupe'
                        print '#' + ' '.join(d1[key].ljust(33) for key in d1.keys())
                        print '#' + ' '.join(d2[key].ljust(33) for key in d2.keys())
                        print '#rm %s' % d1['path']
                        print '#rm %s' % d2['path']
                    else:
                        print '#' + ' '.join(d1[key].ljust(33) for key in d1.keys())
                        print '#' + ' '.join(d2[key].ljust(33) for key in d2.keys())
                        print '#rm %s' % d1['path']
                        print '#rm %s' % d2['path']

                    print

    if options.no_exif:
        for row in c.execute("SELECT * FROM %s WHERE exif_time=? ORDER BY md5" % options.database, ('None',)):
            print '#' + ' '.join(row[key].ljust(33) for key in row.keys())
            #print row
    
    if options.exif_file_mismatch:
        for row in c.execute("SELECT * FROM %s WHERE not exif_time=file_time AND not exif_time=?" % options.database, ('None',)):
            print '#' + ' '.join(row[key].ljust(33) for key in row.keys())
            #print row

    if options.exif_dir_mismatch:
        for row in c.execute("SELECT * FROM %s WHERE not exif_time=dir AND not exif_time=?" % options.database, ('None',)):
            print '#' + ' '.join(row[key].ljust(33) for key in row.keys())
            #print row

    if options.list:
        for row in c.execute("SELECT * FROM %s WHERE path GLOB ?" % options.database, ('*%s*' % options.list,)):
            print '#' + ' '.join(row[key].ljust(33) for key in row.keys())
    

    # Save (commit) the changes
    conn.commit()
    conn.close()
    sys.exit()

else:

    conn = sqlite3.connect('database.sql')
    conn.row_factory = sqlite3.Row
    c = conn.cursor()

    
    if options.action in ['minus', 'rminus', 'xor']:
        if options.action  == 'minus':
            one, two = options.database, options.database2
        else:
            two, one = options.database, options.database2
        c.execute("SELECT md5 FROM %s EXCEPT SELECT md5 from %s" % (one, two))
        res = c.fetchall()
        if not options.action == 'xor':
            sys.stderr.write('subtraction has %d entries\n' % len(res))
            for md5 in res:
                for row in c.execute("SELECT * FROM %s WHERE md5=?" % one, (md5[0],)):
                    print '#' + ' '.join(row[key].ljust(33) for key in row.keys())
        else:
            c.execute("SELECT md5 FROM %s EXCEPT SELECT md5 from %s" % (two, one))
            res.extend(c.fetchall())
            sys.stderr.write('xor has %d entries\n' % len(res))
            for md5 in res:
                for row in c.execute("SELECT * FROM %s WHERE md5=?" % one, (md5[0],)):
                    print '#' + ' '.join(row[key].ljust(33) for key in row.keys())
            
    elif options.action == 'intersection':
        one, two = options.database, options.database2
        c.execute("SELECT md5 FROM %s EXCEPT SELECT md5 from %s" % (one, two))


    exit()


    dbfile = 'database.shelf'
    dbfile2 = 'database2.shelf'

    if os.path.exists(dbfile):
        sys.stderr.write('opening existing database\n')
        #img_dict = shelve.open(dbfile, protocol=-1, writeback=use_writeback)
        #known_paths = shelve.open(pfile, protocol=-1, writeback=use_writeback)
        
        img_shelf1 = shelve.open(dbfile)
        if 'known_paths' in img_shelf1:
            known_paths1 = img_shelf1['known_paths']
        if 'img_dict' in img_shelf1:
            img_dict1 = img_shelf1['img_dict']

    if os.path.exists(dbfile2):
        sys.stderr.write('opening existing database\n')
        #img_dict = shelve.open(dbfile, protocol=-1, writeback=use_writeback)
        #known_paths = shelve.open(pfile, protocol=-1, writeback=use_writeback)
        
        img_shelf2 = shelve.open(dbfile2)
        if 'known_paths' in img_shelf2:
            known_paths2 = img_shelf2['known_paths']
        if 'img_dict' in img_shelf2:
            img_dict2 = img_shelf2['img_dict']

    set1 = set([f for f in flatten_array(img_dict1.values())])
    set2 = set([f for f in flatten_array(img_dict2.values())])

    sys.stderr.write('databases have %d and %d entries\n' % (len(set1), len(set2)))
    
    if options.action == 'intersection':
        newSet = sorted(list(set1 & set2))
        sys.stderr.write('intersection has %d entries\n' % len(newSet))
    elif options.action == 'union':
        newSet = sorted(list(set1 | set2))
        sys.stderr.write('union has %d entries\n' % len(newSet))
    elif options.action == 'xor':
        newSet = sorted(list(set1 ^ set2))
        sys.stderr.write('xor has %d entries\n' % len(newSet))
    elif options.action == 'minus':
        newSet = sorted(list(set1 - set2))
        sys.stderr.write('first - second has %d entries\n' % len(newSet))
    elif options.action == 'rminus':
        newSet = sorted(list(set2 - set1))
        sys.stderr.write('second - first has %d entries\n' % len(newSet))

    for f in newSet:
        if not options.minimal:
            sys.stdout.write('\t%s\n' % f)
        else:
            sys.stdout.write('%s\n' % f.path)


exit()

f = open(sys.argv[1], 'rb')
tags = exifread.process_file(f)

for tag in tags.keys():
    if tag not in ('JPEGThumbnail', 'TIFFThumbnail', 'Filename', 'EXIF MakerNote'):
        print "Key: %s, value %s" % (tag, tags[tag])
exit()
#print '%s' % tags['EXIF DateTimeOriginal']
time = re.sub(':', '-', re.search('[0-9]+:[0-9]+:[0-9]+', '%s' % tags['EXIF DateTimeOriginal']).group(0).strip())
print time


'''
if options.duplicates:
    sys.stderr.write('files with duplicate md5:')
    for md5, group in img_dict.items():
        if len(group) > 1:
            sys.stdout.write('%s\n' % md5)
            for f in group:
                sys.stdout.write('\t%s\n' % f)

if options.exif_file_mismatch:
    sys.stderr.write('files with exif / filedate mismatch:\n')
    #files = []
    files = set()
    for md5, group in img_dict.items():
        for f in group:
            if f.exif_time and f.exif_time != f.file_time:
                #files.append(f)
                files.add(f)
                #sys.stdout.write('\t%s\n' % f)

    if options.reverse:
        files = list(set(img_dict.values()) - files)
    else:
        files = list(files)
        #set([f for f in flatten_array(img_dict1.values())])
    files.sort(key=lambda f:f.path)

    for f in files:
        if not options.minimal:
            sys.stdout.write('\t%s\n' % f)
        else:
            sys.stdout.write('%s\n' % f.path)

if options.exif_dir_mismatch:
    sys.stderr.write('files with exif / directory name mismatch:\n')
    files = []
    for md5, group in img_dict.items():
        for f in group:
            if f.exif_time and f.exif_time != f.dir:
                files.append(f)
                sys.stdout.write('\t%s\n' % f)

    files.sort(key=lambda f:f.path)
    for f in files:
        if not options.minimal:
            sys.stdout.write('\t%s\n' % f)
        else:
            sys.stdout.write('%s\n' % f.path)

if options.no_exif:
    sys.stderr.write('files without exif date:\n')
    files = []
    for md5, group in img_dict.items():
        for f in group:
            if f.exif_time is None:
                files.append(f)

    files.sort(key=lambda f:f.path)
    for f in files:
        if not options.minimal:
            sys.stdout.write('\t%s\n' % f)
        else:
            sys.stdout.write('%s\n' % f.path)

sys.stderr.write('closing database...\n')
if not use_writeback and not options.read_only:
    img_shelf['img_dict'] = img_dict
    img_shelf['known_paths'] = known_paths

'''


