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
            sys.stderr.write('Warning: tree appears to be rooted\n')
    return False

#use argparse module to parse commandline input
parser = argparse.ArgumentParser()

parser.add_argument('-op', '--outgroup-pattern', default=None,
                    help='regex pattern matching taxon label to use as outgroup (single taxon outgroup)')

parser.add_argument('treefiles', nargs='*', default=[], help='nexus treefile to convert')

parser.add_argument('-o', '--outfile', default=None, 
                    help='file to write output to (default is stdout)')

parser.add_argument('-n', '--nexus', action='store_true', default=False, 
                    help='output treefile in nexus rather than newick format')

mut_group = parser.add_mutually_exclusive_group()

mut_group.add_argument('-nb', '--no-bifurcating', action='store_true', default=False, 
                    help='do not include bifurcating trees in output')

mut_group.add_argument('-np', '--no-polytomies', action='store_true', default=False, 
                    help='do not include polytomous trees in output')

mut_group.add_argument('--make-bifurcating', action='store_true', default=False, 
                    help='randomly resolve polytomous nodes with zero-length branches, meaning that all trees will be output and will be bifurcating')

parser.add_argument('--suppress-branchlengths', action='store_true', default=False, 
                    help='strip branchlengths from output trees (default False)')

parser.add_argument('-p', '--prune-patterns', action='append', default=None, 
                    help='patterns for taxon names to strip from trees before output.  Single pattern per flag, but can appear multiple times')

options = parser.parse_args()

#sys.stderr.write('Reading %s ...\n' % options.treefiles)

intrees = dendropy.TreeList()
if not options.treefiles:
    sys.stderr.write('NOTE: reading trees from stdin')
    trees = sys.stdin.read()
    try:
        intrees.extend(dendropy.TreeList.get_from_string(trees, "nexus"))
    except dendropy.error.DataError:
        intrees.extend(dendropy.TreeList.get_from_string(trees, "newick"))

else:
    for tf in options.treefiles:
        #try two input formats
        try:
            intrees.extend(dendropy.TreeList.get_from_path(tf, "nexus"))
        except dendropy.error.DataError:
            intrees.extend(dendropy.TreeList.get_from_path(tf, "newick"))

sys.stderr.write('read %d trees\n' % len(intrees))

#treestr = '(O._barthii_AA:0.00157155,(((O._brachyantha_FF:0.10458481,(O._punctata_BB:0.00266559,O._minuta_BB:0.01210456):0.01556435):0.00268608,(O._officinalis_CC:0.10078888,O._minuta_CC:0.02347313):0.01668656):0.03394209,((O._sativaj_AA:0.01511099,O._rufipogon_AA:0.00251092):0.00401496,O._nivara_AA:0.002933):0.00296048):0.00068407,O._glaberrima_AA:1e-08);'
#intree = dendropy.Tree()
#intree.read_from_string(treestr, 'newick')

out = open(options.outfile, 'w') if options.outfile else sys.stdout
log = sys.stderr

outtrees = dendropy.TreeList()
ignoredCount = 0
outgroupIgnoredCount = 0
madeBifurcating = 0
for intree in intrees:
    hasPoly = check_for_polytomies(intree)
    if options.no_bifurcating and not hasPoly:
        ignoredCount += 1
        #log.write('ignoring bifurcating tree\n')
    elif options.no_polytomies and hasPoly:
        ignoredCount += 1
        #log.write('ignoring polytomous tree\n')
    else:
        if options.make_bifurcating and hasPoly:
            intree.resolve_polytomies(update_splits=True)
            madeBifurcating += 1
        to_remove = set()
        #prune taxa first with patterns, THEN look for an outgroup pattern.
        #outgroup pattern could be specified that matches something that has
        #already been deleted
        if options.prune_patterns is not None:
            leaves = intree.leaf_nodes()
            for l in leaves:
                for to_prune in options.prune_patterns:
                    if re.search(to_prune, l.taxon.label) is not None:
                        to_remove.add(l.taxon.label)
                        break
            to_retain = set(l.taxon.label for l in leaves) - to_remove
            #intree.retain_taxa_with_labels(labels=to_retain)
            intree.prune_taxa_with_labels(labels=to_remove)
            #these are called on TreeLists - not sure if applicable here
            intree.taxon_set = intree.infer_taxa()
            intree.reindex_subcomponent_taxa()

        if options.outgroup_pattern is not None:
            outgroup = None
            leaves = intree.leaf_nodes()
            for l in leaves:
                if re.search(options.outgroup_pattern, l.taxon.label) is not None:
                    outgroup = l
                    break

            if outgroup is None:
                #log.write('ignoring tree without specified outgroup\n')
                outgroupIgnoredCount += 1
                continue
            else:
                #if the tree was already rooted, this will remove that root node
                #rooting halves the branchlength of the chosen branch
                if outgroup.edge_length:
                    intree.reroot_at_edge(outgroup.edge, length1=outgroup.edge_length / 2.0, length2=outgroup.edge_length / 2.0, update_splits=False, delete_outdegree_one=True) 
                else:
                    intree.reroot_at_edge(outgroup.edge, update_splits=False, delete_outdegree_one=True) 
        
        outtrees.append(intree)

if ignoredCount > 0:
    log.write('ignored %d trees\n' % ignoredCount)
if outgroupIgnoredCount > 0:
    log.write('ignored %d trees because of missing outgroup\n' % outgroupIgnoredCount)
if madeBifurcating > 0:
    log.write('%d polytomous trees arbitrarily resolved\n' % madeBifurcating)

if outtrees:
    log.write('writing %d trees\n' % len(outtrees))
    if options.nexus:
        outtrees.write(out, "nexus", suppress_edge_lengths=options.suppress_branchlengths)
    else:
        outtrees.write(out, "newick", suppress_edge_lengths=options.suppress_branchlengths)
else:
    log.write('no trees to output?\n')

