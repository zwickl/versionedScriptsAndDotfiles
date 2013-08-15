#!/usr/bin/env python

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
            #if item.is_rooted != t.is_rooted or numNodes == len(t.nodes()):
                #print 'testing'
            if treecalc.symmetric_difference(t, item) == 0:
                #print 'same'
                return True
            '''
            print 'testing'
            if treecalc.symmetric_difference(t, item) == 0:
                print 'same'
                return True
            '''
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

    def masked_frequency_of_splitlist(self, **kwargs):
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
        requiredSplits = len(targetSplits)
        for tree in self:
            '''
            #print dir(tree)
            #print tree.as_newick_string()
            for internal in tree.postorder_internal_node_iter():
                print dir(internal)
                print internal.label
            for edge in tree.postorder_edge_iter():
                print dir(tree.seed_node)
                print tree.seed_node.description()
            '''
            if not hasattr(tree, "split_edges"):
                treesplit.encode_splits(tree)
            total += 1
            matchedSplits = 0
            incompatible = False
            #work through the required splits
            for targetSplit in targetSplits:
                compSplit = (~targetSplit & partialMask)
                #work through the splits in this tree
                for test_split in tree.split_edges:
                    #mask out unimportant taxa
                    masked_test = (test_split & partialMask)
                    if not treesplit.is_compatible(test_split, targetSplit, partialMask):
                        incompatible = True
                        break
                    if targetSplit == masked_test or compSplit == masked_test :
                        matchedSplits += 1
                        break
                if incompatible:
                    break
            if not incompatible and matchedSplits == requiredSplits:
                found += 1
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


