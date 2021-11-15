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
            #data = fh.read(8192)
            data = fh.read(16384)
            if not data:
                break
            m.update(data)
    return m.hexdigest()


class FileInfo(object):
    def __init__(self, fname):
        #get rid of multiple slashes in path!
        fname = re.sub('/+', '/', fname)
        self.path = fname
        self.fname = os.path.basename(fname)
        self.dir = fname.split('/')[-2]
        #self.stat = os.stat(fname)
        #self.file_time = time.ctime(os.path.getctime(fname))
        self.file_time = str(date.fromtimestamp(os.path.getmtime(fname)))
        self.md5 = md5Checksum(fname)
        self.hash = hash(self.md5)
        self.size = os.path.getsize(fname)

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
        self.altname = re.sub('mov$', 'thm', re.sub('MOV$', 'THM', fname))
        if not os.path.isfile(self.altname):
            self.altname = None
        with open(self.altname or fname, 'rb') as f:
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
            try:
                self.epoch_seconds = time.mktime(time.strptime(self.tags['EXIF DateTimeOriginal'].values, '%Y:%m:%d %H:%M:%S'))
            except ValueError:
                try:
                    self.epoch_seconds = time.mktime(time.strptime(self.tags['EXIF DateTimeOriginal'].values, '%Y:%m:%d %H:%M: %S'))
                except ValueError:
                    sys.stderr.write('Could not parse exif data %s: %s\n' % (fname, self.tags['EXIF DateTimeOriginal']))
                    self.epoch_seconds = None

    def set_timestamp_to_exif(self):
        os.utime(self.path, (self.epoch_seconds, self.epoch_seconds))

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
        return (self.path, self.dir, self.fname, self.file_time, ex_time, self.size, self.md5)

    def header(self):
        return 'path\texif_time\tfile_time\tdir'


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


def make_tabular_string(img_dict):
    print '#' + ' '.join(str(img_dict[key]).ljust(width) for key, width in zip(img_dict.keys(), [72, 20, 26, 12, 12, 12])), match_string(img_dict)


def make_tabular_header(img_dict):
    print '#' + ' '.join(key.ljust(width) for key, width in zip(img_dict.keys(), [72, 20, 26, 12, 12, 12]))


def element_identical(first, sec, el):
    return first[el] == sec[el]


def identical_elements_string(first, sec):
    ident = []
    for el in ['path', 'dir', 'filename', 'file_time', 'exif_time', 'md5', 'size']:
        if element_identical(first, sec, el):
            ident.append(el)
    return ' '.join(ident)


def match_string(entry):
    if entry['exif_time'] == entry['file_time'] == entry['dir']:
        return 'EFD'
    elif entry['exif_time'] == entry['file_time'] == re.sub('_', '-', entry['dir']):
        return 'EFd'
    elif entry['exif_time'] == entry['file_time']:
        return 'EF'
    elif entry['file_time'] == entry['dir']:
        return 'FD'
    elif entry['file_time'] == re.sub('_', '-', entry['dir']):
        return 'Fd'
    elif entry['exif_time'] == entry['dir']: 
        return 'ED'
    elif entry['exif_time'] == re.sub('_', '-', entry['dir']):
        return 'Ed'
    return '-'


def table_report(cursor):
    cursor.execute('SELECT * FROM sqlite_master WHERE type = "table"')
    tab_names = [ tab[1] for tab in cursor.fetchall() ]
    for tab in tab_names:
        cursor.execute('SELECT COUNT(*) AS num FROM %s' % tab)
        total = cursor.fetchone()[0]
        print tab, total


parser = argparse.ArgumentParser(description='Create a database and manipulate image files')

parser.add_argument('-d', '--duplicates', action='store_true', default=False, help='output list of duplicate files to stdout, with bash script to remove them')
parser.add_argument('--move-dir', default=None, help='path to move duplicate files to in bash script, rather than delete them')

parser.add_argument('--duplicate-key', default='md5', help='choose the field used to determine duplicates (default md5, also filename, size, file_time, exif_time)')

parser.add_argument('-v', '--reverse', action='store_true', default=False, help='reverse sense of no-exif or mismatch match (files NOT meeting condition)')

parser.add_argument('--no-exif', action='store_true', default=False, help='output list of files without exif information to stdout')

parser.add_argument('--exif-file-mismatch', action='store_true', default=False, help='output list of files with conflict between file date and exif info to stdout')

parser.add_argument('--exif-dir-mismatch', action='store_true', default=False, help='output list of files with conflict between directory name and exif info to stdout')

parser.add_argument('--match-file-to-exif', action='store_true', default=False, help='change files modify and access times to match exif create time')

parser.add_argument('--read-only', action='store_true', default=False, help='perform actions on the database (even adding files), but don\'t write it back to file (faster)')

parser.add_argument('--minimal', action='store_true', default=False, help='only output a list of filenames for the chosen option, no other information')

parser.add_argument('-db', '--database', default='database.sql', help='name of database file to load/create (default database.sql)')

parser.add_argument('-t', '--table', default='table1', help='name of table in database to load/create (default is "table1")')
parser.add_argument('-t2', '--table2', default=None, help='name of 2nd table in database use in table comparisons with --action')

parser.add_argument('--action', default=None, choices=['intersection', 'union', 'xor', 'minus', 'rminus'], help='set operations to perform on two databases')

parser.add_argument('--list-tables', action='store_true', default=False, help='list names of tables currently in database')
parser.add_argument('--copy-table', action='store_true', default=False, help='copy contents of table into table2, creating table2 if necessary')
parser.add_argument('--drop-table', action='store_true', default=False, help='delete the table specified with --table, or table1 by default')

parser.add_argument('-l', '--list', action='store_true', default=False, help='list entries in table, possibly matching --pattern')

parser.add_argument('--pattern', default="*", help='list entries containing this pattern (NOT regex, but *''s and ?\'s needed to match entire path)')

parser.add_argument('files', nargs='*', default=[], help='image or movie files to process')

options = parser.parse_args()

use_writeback = False

if options.list_tables:
    conn = sqlite3.connect(options.database)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()
    table_report(c)

elif options.drop_table:
    conn = sqlite3.connect(options.database)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()
    c.execute('DROP TABLE %s' % options.table)
    print 'remaining tables:'
    table_report(c)

elif options.copy_table:
    conn = sqlite3.connect(options.database)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()
    c.execute('''CREATE TABLE if not exists %s
            (path TEXT PRIMARY KEY, dir TEXT, filename TEXT, file_time TEXT, exif_time TEXT, size INTEGER,  md5 TEXT)''' % options.table2)
    c.execute('''SELECT * from %s''' % (options.table))
    entries = c.fetchall()
    c.executemany('''INSERT OR IGNORE INTO %s VALUES (?, ?, ?, ?, ?, ?, ?)''' % options.table2, [ tuple(ins) for ins in entries ])
    # Save (commit) the changes
    conn.commit()
    conn.close()
    sys.exit()

elif not options.action:
    #img_shelf, img_dict, known_paths = open_image_shelf(options.database)

    conn = sqlite3.connect(options.database)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE if not exists %s
                (path TEXT PRIMARY KEY, dir TEXT, filename TEXT, file_time TEXT, exif_time TEXT, size INTEGER,  md5 TEXT)''' % options.table)
    except sqlite3.OperationalError:
        print 'Problem creating/accessing table'
        print options.table


    img_dict = {}
    
    numAdded = 0
    numHad = 0
    for fn in options.files:
        if os.path.isfile(fn):
                
            c.execute('SELECT * FROM %s WHERE path=?' % options.table, (fn,))
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

    c.executemany("INSERT OR IGNORE INTO %s VALUES (?, ?, ?, ?, ?, ?, ?)" % options.table, [ img.tuple_repr() for img in flatten_array(img_dict.values()) ])
 
    c.execute('SELECT COUNT(*) AS num FROM %s' % options.table)
    total = c.fetchone()[0]
    sys.stderr.write('\n%d items added to database\n%d items were already in database\n%d total items in database\n' % (numAdded, numHad, total))
    
    if options.duplicates:
        field = options.duplicate_key
        #c.execute("SELECT md5 FROM %s WHERE path GLOB ? GROUP BY md5 HAVING count(md5) > 1" % options.table, (options.pattern,))
        c.execute("SELECT %s FROM %s WHERE path GLOB ? GROUP BY %s HAVING count(%s) > 1" % (field, options.table, field, field), (options.pattern,))
        first = True
        for vals in c.fetchall():
            #for row in c.execute("SELECT * FROM %s WHERE md5=?" % options.database, (md5[0],)):
            #    print row
            #print
            #dupes = [ row for row in c.execute("SELECT * FROM %s WHERE md5=?" % options.table, (md5[0],)) ]
            dupes = [ row for row in c.execute("SELECT * FROM %s WHERE %s=?" % (options.table, field), (vals[0],)) ]
            #print '\n'.join([str(d) for d in dupes])
            for n1, d1 in enumerate(dupes):
                for d2 in dupes[n1 + 1:]:
                    if first:
                        make_tabular_header(d1)
                        first = False
                    if list(d1)[1:] == list(d2)[1:]:
                        print '#all but path dupe'
                    elif list(d1)[2:] == list(d2)[2:]:
                        print '#all but directory dupe'
                    elif list(d1)[3:] == list(d2)[3:]:
                        print '#exif-md5 dupe'
                    sys.stdout.write('#%s\n' % identical_elements_string(d1, d2))
                    make_tabular_string(d1)
                    make_tabular_string(d2)

                    plen1 = len(d1['path']) 
                    plen2 = len(d2['path'])
                    flen1 = len(d1['filename']) 
                    flen2 = len(d2['filename'])
                    
                    #first_shorter = len1 < len2
                    #first_alpha = d1['path'] < d2['path']

                    #first favor pics where everything matches
                    remove_first = None
                    if (len(match_string(d1)) == 3 and 'D' in match_string(d1)) and (len(match_string(d2)) == 3 and 'D' in match_string(d2)) and len1 != len2:
                        remove_first = False if first_shorter else True
                    elif len(match_string(d1)) == 3 and 'D' in match_string(d1):
                        remove_first = False
                    elif len(match_string(d2)) == 3 and 'D' in match_string(d2):
                        remove_first = True
                    #lowercase d means that the directory matches except for _ vs - differences
                    elif len(match_string(d1)) == 3 and 'd' in match_string(d1):
                        remove_first = False
                    elif len(match_string(d2)) == 3 and 'd' in match_string(d2):
                        remove_first = True
                    elif d1['exif_time'] not in [None, 'None'] and d2['exif_time'] in [None, 'None']:
                        remove_first = False
                    elif d1['exif_time'] in [None, 'None'] and d2['exif_time'] not in [None, 'None']:
                        remove_first = True
                    elif d1['exif_time'] not in [None, 'None'] and d1['exif_time'] != d1['dir']:
                        remove_first = True
                    elif d2['exif_time'] not in [None, 'None'] and d2['exif_time'] != d2['dir']:
                        remove_first = False
                    elif 'temp' in d1['path'].lower():
                        remove_first = True
                    elif 'temp' in d2['path'].lower():
                        remove_first = False
                    else:
                        use_fnamelen = True
                        
                        #then favor the shorter or longer path or filename
                        if use_fnamelen:
                            if flen1 != flen2:
                                if flen1 > flen2:
                                    remove_first = False
                                else:
                                    remove_first = True
                        else:
                            if plen1 != plen2:
                                if plen1 < plen2:
                                    remove_first = False
                                else:
                                    remove_first = True
                    
                    '''
                    elif d1['path'] < d2['path']:
                        #finally, favor the alphanumerically first
                        remove_first = False
                    else:
                        remove_first = True
                    '''

                    if options.move_dir:
                        target1 = os.path.join(options.move_dir, d1['dir'])
                        target2 = os.path.join(options.move_dir, d2['dir'])
                        if remove_first:
                            print 'mkdir -p \'%s\'; mv -n \'%s\' \'%s\'/' % (target1, d1['path'], target1)
                            print '#mkdir -p \'%s\'; mv -n \'%s\' \'%s\'/' % (target2, d2['path'], target2)
                        else:
                            print '#mkdir -p \'%s\'; mv -n \'%s\' \'%s\'/' % (target1, d1['path'], target1)
                            print 'mkdir -p \'%s\'; mv -n \'%s\' \'%s\'/' % (target2, d2['path'], target2)
                    else:
                        if remove_first:
                            print 'rm \'%s\'' % d1['path']
                            print '#rm \'%s\'' % d2['path']
                        else:
                            print '#rm \'%s\'' % d1['path']
                            print 'rm \'%s\'' % d2['path']
                    print

    else:
        first = True
        if options.no_exif:
            for row in c.execute("SELECT * FROM %s WHERE path GLOB ? AND exif_time=? ORDER BY md5" % options.table, (options.pattern, 'None',)):
                if first:
                    make_tabular_header(row)
                    first = False
                make_tabular_string(row)
        
        if options.exif_file_mismatch or options.match_file_to_exif:
            #for row in c.execute("SELECT * FROM %s WHERE not exif_time=file_time AND not exif_time=?" % options.table, ('None',)): 
            rows = sorted([ r for r in c.execute("SELECT * FROM %s WHERE path GLOB ? AND not exif_time=file_time AND not exif_time=?" % options.table, (options.pattern, 'None',)) if r['file_time'] > r['exif_time'] and r['exif_time'] == r['dir'] ], key=lambda el: (el['exif_time'], el['filename']))
            yn = None
            for row in rows:
                if first:
                    make_tabular_header(row)
                    first = False
                make_tabular_string(row)
                if options.match_file_to_exif:
                    if yn != 'A':
                        yn = raw_input('Change file timestamp to match exif?')
                    if yn.lower() == 'y' or yn.lower() == 'a':
                        f = ImageFile(row['path'])
                        f.set_timestamp_to_exif()
                        newf = ImageFile(f.path)
                        print newf.file_time
                        c.execute("UPDATE %s SET file_time=? WHERE path=?" % options.table, (newf.file_time, newf.path))

        if options.exif_dir_mismatch:
            for row in c.execute("SELECT * FROM %s WHERE path GLOB ? AND not exif_time=dir AND not exif_time=?" % options.table, (options.pattern, 'None',)):
                if first:
                    make_tabular_header(row)
                    first = False
                make_tabular_string(row)

        if options.list:
            for row in c.execute("SELECT * FROM %s WHERE path GLOB ?" % options.table, (options.pattern,)):
                if first:
                    make_tabular_header(row)
                    first = False
                make_tabular_string(row)

    # Save (commit) the changes
    conn.commit()
    conn.close()
    sys.exit()

else:
    conn = sqlite3.connect(options.database)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()
    
    if options.action in ['minus', 'rminus', 'xor']:
        if options.action  == 'minus':
            one, two = options.table, options.table2
        else:
            two, one = options.table, options.table2
        c.execute("SELECT %s FROM %s EXCEPT SELECT %s from %s" % (options.duplicate_key, one, options.duplicate_key, two))
        res = c.fetchall()
        if not options.action == 'xor':
            sys.stderr.write('subtraction has %d entries\n' % len(res))
            for val in res:
                for row in c.execute("SELECT * FROM %s WHERE %s=?" % (one, options.duplicate_key), (val[0],)):
                    make_tabular_string(row)
        else:
            c.execute("SELECT md5 FROM %s EXCEPT SELECT md5 from %s" % (two, one))
            res.extend(c.fetchall())
            sys.stderr.write('xor has %d entries\n' % len(res))
            for md5 in res:
                for row in c.execute("SELECT * FROM %s WHERE md5=?" % one, (md5[0],)):
                    make_tabular_string(row)
            
    elif options.action == 'intersection':
        one, two = options.table, options.table2
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


