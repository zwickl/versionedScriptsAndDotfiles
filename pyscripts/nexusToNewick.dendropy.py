#!/usr/bin/env python

import sys
import re
import argparse
import dendropy

def check_for_polytomies(tree):
    for node in tree.postorder_node_iter():
        #sys.stderr.write("%s\n" % len(node.adjacent_nodes()))
        if len(node.adjacent_nodes()) > 3:
            return True
        elif len(node.adjacent_nodes()) == 2:
            sys.stderr.write('Warning: tree appears to be rooted')
    return False

#use argparse module to parse commandline input
parser = argparse.ArgumentParser()

parser.add_argument('-op', '-outgroup-pattern', dest='patt', required=False, default=None,
                    help='regex pattern matching taxon label to use as outgroup (single taxon outgroup)')

parser.add_argument('treefiles', nargs='*', default=[], help='nexus treefile to convert')

parser.add_argument('-o', '--outfile', dest='outfile', required=False, default=None, 
                    help='file to write output to (default is stdout)')

parser.add_argument('-n', '--nexus', dest='outputNexus', required=False, default=False, 
                    help='output treefile in nexus rather than newick format')

parser.add_argument('-nb', '--no-bifurcating', dest='ignoreBif', action='store_true', required=False, default=False, 
                    help='do not include bifurcating trees in output')

parser.add_argument('-np', '--no-polytomies', dest='ignorePoly', action='store_true', required=False, default=False, 
                    help='do not include polytomous trees in output')


parsed = parser.parse_args()

#sys.stderr.write('Reading %s ...\n' % parsed.treefiles)

intrees = dendropy.TreeList()
for tf in parsed.treefiles:
    intrees.append(dendropy.Tree.get_from_path(tf, "nexus"))
sys.stderr.write('read %d trees\n' % len(intrees))

#treestr = '(O._barthii_AA:0.00157155,(((O._brachyantha_FF:0.10458481,(O._punctata_BB:0.00266559,O._minuta_BB:0.01210456):0.01556435):0.00268608,(O._officinalis_CC:0.10078888,O._minuta_CC:0.02347313):0.01668656):0.03394209,((O._sativaj_AA:0.01511099,O._rufipogon_AA:0.00251092):0.00401496,O._nivara_AA:0.002933):0.00296048):0.00068407,O._glaberrima_AA:1e-08);'
#intree = dendropy.Tree()
#intree.read_from_string(treestr, 'newick')

if parsed.outfile:
    out = open(parsed.outfile, 'w')
else:
    out = sys.stdout

outtrees = dendropy.TreeList()
ignoredCount = 0
outgroupIgnoredCount = 0
for intree in intrees:
    hasPoly = check_for_polytomies(intree)
    if parsed.ignoreBif is True and hasPoly is False:
        ignoredCount += 1
        #sys.stderr.write('ignoring bifurcating tree\n')
    elif parsed.ignorePoly is True and hasPoly is True:
        ignoredCount += 1
        #sys.stderr.write('ignoring polytomous tree\n')
    else:
        if parsed.patt is not None:
            outgroup = None
            leaves = intree.leaf_nodes()
            for l in leaves:
                if re.search(parsed.patt, l.taxon.label) is not None:
                    outgroup = l
                    break

            if outgroup is None:
                #sys.stderr.write('ignoring tree without specified outgroup\n')
                outgroupIgnoredCount += 1
                continue
            else:
                #if the tree was already rooted, this will remove that root node
                #rooting halves the branchlength of the chosen branch
                intree.reroot_at_edge(outgroup.edge, length1=outgroup.edge_length / 2.0, length2=outgroup.edge_length / 2.0, update_splits=False, delete_outdegree_one=True) 
        outtrees.append(intree)

if ignoredCount > 0:
    sys.stderr.write('ignored %d trees\n' % ignoredCount)
if outgroupIgnoredCount > 0:
    sys.stderr.write('ignored %d trees because of missing outgroup\n' % outgroupIgnoredCount)
if len(outtrees) > 0:
    sys.stderr.write('writing %d trees\n' % len(outtrees))
    if parsed.outputNexus:
        outtrees.write(out, "nexus")
    else:
        outtrees.write(out, "newick")
else:
    sys.stderr.write('no trees to output?')

