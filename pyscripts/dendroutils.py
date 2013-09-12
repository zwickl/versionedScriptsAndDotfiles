#!/usr/bin/env python

import copy
from dendropy import Tree, TreeList, treesplit

class MyTree(Tree):
    '''This is my override of denodropy Tree, which redefines equality as 
    '''
    def __eq__(self, other):
        return treecalc.symmetric_difference(self, other) == 0
    def __hash__(self):
        h = hash(''.join(sorted(self.as_newick_string())))
        print h
        return h

class MyTreeList(TreeList):

    def __contains__(self, item):
        '''overridden function to allow basic use of 'in' keyword, like
        if treeX in treelistY: 
            blah
        NOT very efficient
        '''
        numNodes = len(item.nodes())
        for t in self:
            #don't do the expensive symmetric dist if # of nodes isn't same
            #although safest to do it if the trees differ in their rootednes

            #EDIT - this is too dangerous - dendropy can consider trees unrooted
            #even if they have a bifurcating root.  Just test.
            if treecalc.symmetric_difference(t, item) == 0:
                return True
        #print 'no match'
        return False

    def frequency_of_identical_trees(self, targetTree):
        numNodes = len(targetTree.nodes())
        count = 0
        total = len(self)

        for tree in self:
            #don't do the expensive symmetric dist if # of nodes isn't same
            #although safest to do it if the trees differ in their rootedness
            #EDIT - this is too dangerous - dendropy can consider trees unrooted
            #even if they have a bifurcating root.  Just test.
            #if targetTree.is_rooted != t.is_rooted or numNodes == len(t.nodes()):
            if treecalc.symmetric_difference(tree, targetTree) == 0:
                count += 1

        return float(count) / total

    def masked_frequency_of_split(self, **kwargs):
        """
        DJZ - this is my adaptation of frequency_of_split that takes a
        taxon mask, in my own derived TreeList class

        Given a split or bipartition specified as:

            - a split bitmask given the keyword 'split_bitmask'
            - a list of `Taxon` objects given with the keyword `taxa`
            - a list of taxon labels given with the keyword `labels`
            - a list of oids given with the keyword `oids`

        this function returns the proportion of trees in self
        in which the split is found.
        """
#DJZ
        partialMask = kwargs["mask"] if "mask" in kwargs else self.taxon_set.all_taxa_bitmask()

        if "split_bitmask" in kwargs:
            targetSplit = kwargs["split_bitmask"]
        else:
            targetSplit = self.taxon_set.get_taxa_bitmask(**kwargs)
            k = kwargs.values()[0]
            if treesplit.count_bits(targetSplit) != len(k):
                raise IndexError('Not all taxa could be mapped to split (%s): %s' \
                    % (self.taxon_set.split_bitmask_string(targetSplit), k))
        found = 0
        total = 0
        for tree in self:
            if not hasattr(tree, "split_edges"):
                treesplit.encode_splits(tree)
            total += 1
            #if split in tree.split_edges:
            #    found += 1

            compSplit = (~targetSplit & partialMask)
            for test_split in tree.split_edges:
                if not treesplit.is_compatible(test_split, targetSplit, partialMask):
                    break
                masked_test = (test_split & partialMask)
                if targetSplit == masked_test or compSplit == masked_test :
                    found += 1
                    break

        return float(found)/total

    def masked_frequency_of_splitlist(self, returnMatches=False, **kwargs):
        """
        DJZ - this is my adaptation of frequency_of_masked_split that takes a
        taxon mask, in my own derived TreeList class

        Given a LIST of splits or bipartitions specified as:

            - a split bitmask given the keyword 'split_bitmask'
            - a list of `Taxon` objects given with the keyword `taxa`
            - a list of taxon labels given with the keyword `labels`
            - a list of oids given with the keyword `oids`

        this function returns the proportion of trees in self
        in which all of the splits are found.
        NOTE: This is not that useful in some cases here you call it sucessively with
        different numbers of splits and expect the freqs to add up to 1.0
        """
        
        #if returnMatches is requested, return matching trees
        matches = []

        partialMask = kwargs["mask"] if "mask" in kwargs else self.taxon_set.all_taxa_bitmask()

        if "split_bitmask" in kwargs:
            targetSplits = kwargs["split_bitmask"]
        else:
            split = self.taxon_set.get_taxa_bitmask(**kwargs)
            k = kwargs.values()[0]
            if treesplit.count_bits(split) != len(k):
                raise IndexError('Not all taxa could be mapped to split (%s): %s' \
                    % (self.taxon_set.split_bitmask_string(split), k))

        found = 0
        total = 0
        for tnum, tree in enumerate(self):
            if not hasattr(tree, "split_edges"):
                treesplit.encode_splits(tree)
            total += 1
            matchedSplits = 0
            incompatible = False
            #work through the required splits
            for num, targetSplit in enumerate(targetSplits):
                compSplit = (~targetSplit & partialMask)
                #work through the splits in this tree
                for test_split in tree.split_edges:
                    #mask out unimportant taxa
                    masked_test = (test_split & partialMask)
                    #don't need to test anything if masked_test is empty (i.e., no taxa in partialMask appear on opposite sides
                    #of test_split
                    if masked_test:
                        #print '%13s %13s %13s %d' % (bin(targetSplit), bin(test_split), bin(masked_test), masked_test),
                        if not treesplit.is_compatible(test_split, targetSplit, partialMask):
                            incompatible = True
                            break
                        if targetSplit == masked_test or compSplit == masked_test:
                            matchedSplits += 1
                            break
                if incompatible:
                    break
            if not incompatible and matchedSplits == len(targetSplits):
                found += 1
                if returnMatches:
                    matches.append(copy.deepcopy(tree))
        if returnMatches:
            return float(found)/total, matches
        else:
            return float(found)/total

    def generate_all_trees_for_taxon_list(self, taxon_list, min_bipartitions=None, max_bipartitions=None, criterion=None):
        '''Will call functions to generate newick strings representing all possible trees for taxon set.  Work must be done here to make that list
        fully unique allowing for differences in tree rotation.  Can pass min and max bipartitions to control resolvedness of trees, or omit to 
        only generate fully resolved
        This is impractically slow for > 6 taxa.
        >>> TL = MyTreeList()
        >>> TL.generate_all_trees_for_taxon_list(['a', 'b', 'c', 'd'])
        >>> len(TL)
        3
        >>> TL = MyTreeList()
        >>> TL.generate_all_trees_for_taxon_list(['a', 'b', 'c', 'd', 'e'])
        >>> len(TL)
        15

        #don't think that this works for trees with polytomies
        #>>> TL = MyTreeList()
        #>>> TL.generate_all_trees_for_taxon_list(['a', 'b', 'c', 'd', 'e'], min_bipartitions=0, max_bipartitions=2)
        #>>> len(TL)
        #26
        
        #>>> TL = MyTreeList()
        #>>> TL.generate_all_trees_for_taxon_list(['a', 'b', 'c', 'd', 'e', 'f'])
        #>>> len(TL)
        105
        #>>> TL = MyTreeList()
        #>>> TL.generate_all_trees_for_taxon_list(['a', 'b', 'c', 'd', 'e', 'f', 'g'])
        #>>> len(TL)
        945
        #>>> TL = MyTreeList()
        #>>> TL.generate_all_trees_for_taxon_list(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'])
        #>>> len(TL)
        10395
        '''
        ntax = len(taxon_list)

        #default to fully resolved trees
        min_bipartitions = max(min_bipartitions,  0) if min_bipartitions else ntax - 3
        max_bipartitions = max(max_bipartitions,  ntax - 3) if max_bipartitions else ntax - 3

        componentLists = combine_components_and_uniqueify(taxon_list, min_components=ntax - max_bipartitions, max_components=ntax - min_bipartitions, criterion=criterion)

        fullString = ''
        for componentList in componentLists:
            #the component list is a list of tuples, so actually the repr is exactly newick representation
            fullString += repr(componentList) + ';'
        self.read_from_string(fullString, 'newick')

        #print len(self)
        #for t in self:
        #    print t.as_newick_string()
        newList = TreeList(self, taxon_set=self.taxon_set)
        #print len(newList)
        self[:] = []
        for num, tr in enumerate(newList):
            #if tr not in TreeList(self[num+1:]):
            if tr not in self:
                self.append(tr)
                #print tr.as_newick_string()


#this was a hack to ensure that only single taxa were combined, using combine_components, 
#which works around multiple represenations of same tree, but only for 4 or 5 taxa
#this has been deprecated
def arguments_not_list_or_tuple(one, two):
    for t in [list, tuple]:
        if isinstance(one, t) or isinstance(two, t):
            return False
    return True


def combine_components(seedComp, min_components=1, max_components=20, criterion=None):
    '''This will recursively combine components to make all possible nested tuples of orginal list
    Written generally, but not clear what the use would be besides generating newick strings to represent
    trees.
    NOTE: not all trees are really unique due to equality of differently oriented trees.  Need to filter
    returned list to make unique.
    '''
    seedLevel = len(seedComp) 
    if seedLevel < min_components:        return []
    elif seedLevel == min_components:
        return [seedComp]

    toReturn = []    
    if seedLevel <= max_components:
        toReturn.append(seedComp)    
        for num1 in xrange(seedLevel - 1):
            for num2 in xrange(num1 + 1, seedLevel):
                if not criterion or criterion(seedComp[num1], seedComp[num2]):
                    #print seedComp[num1], seedComp[num2]
                    newList = list(seedComp)
                    newItem = [newList.pop(num2), newList.pop(num1)]
                    newList.append(tuple(sorted(newItem)))
                    newList.sort()
                    toReturn.extend(combine_components(newList, min_components=min_components, max_components=max_components, criterion=criterion))
        
    return toReturn


def combine_components_and_uniqueify(seedComponents, min_components=3, max_components=20, criterion=None):
    '''Intent here was to filter what is returned from combine_components such that all represented 
    unique trees.  It does not entirely do so.  It does NOT completely uniquify things if trees 
    are oriented differently or nodes rotated. Use of a criterion like arguments_not_list_or_tuple 
    was a previous hack to allow generation of unique trees in the special case of < 6 taxa.  Not needed now.
    '''
    compLists = combine_components(seedComponents, min_components=min_components, max_components=max_components, criterion=criterion)
    compTuples = [ tuple(c) for c in compLists ]
    compSet = set(compTuples)
    compSet = list(compSet)
    #sort by resolvedness
    compSet.sort(key=lambda t:len(t), reverse=True)

    return compSet


def generate_quartet_lists_from_file(filename='quartetList'):
    '''this generates the possible quartets of taxa from the formatted taxon strings, like
    one, three, four, five
    two, three, four, five
    Or, you can use colons to indicate groups of taxa to make quartets combinatorically,
    i.e. this would be equivalent to the two lines above:
    one:two, three, four, five
    '''

    allQuartets = set()
    comboQuartets = []

    for line in open(filename, 'rb'):
        #ignore blank lines and "comments" starting with #
        if len(line.strip()) > 0 and line[0] != '#':
            if ':' in line:
                comboQuartets.append([tax.strip() for tax in line.split(',')])
            else:
                allQuartets |= set([tuple([tax.strip() for tax in line.split(',')])])

    quartsFromCombos = []
    for combo in comboQuartets:
        thisCombo = []
        for subset in combo:
            thisSubset = []
            for el in subset.split(':'):
                if thisCombo:
                    thisSubset.extend([ p + [el] for p in thisCombo ])
                else:
                    thisSubset.append([el])
            thisCombo = thisSubset
        quartsFromCombos.extend(thisCombo)
    allQuartets |= set([tuple(q) for q in quartsFromCombos])

    #allQuartets was initially a set of tuples, to ensure that duplicate quarts weren't kept, now
    #convert back to list of lists
    return [list(q) for q in allQuartets]


