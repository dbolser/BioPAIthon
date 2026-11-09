# Copyright (C) 2013 by Yanbo Ye (yeyanbo289@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the Bio.Phylo.Consensus module."""

import os
import random
import tempfile
import unittest
from io import StringIO

import support

from Bio import Align
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo import BaseTree
from Bio.Phylo import Consensus
from Bio.Phylo.Consensus import _BitString
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

temp_dir = tempfile.mkdtemp()


def read_newick(text):
    return Phylo.read(StringIO(text), "newick")


def inner_clades(tree):
    """Map the tip names of each clade below the root to its confidence."""
    return {
        frozenset(term.name for term in clade.get_terminals()): clade.confidence
        for clade in tree.find_clades(terminal=False)
        if clade is not tree.root
    }


class BitStringTest(unittest.TestCase):
    """Test for _BitString class."""

    def test_bitstring(self):
        bitstr1 = _BitString("0011")
        bitstr2 = _BitString("0101")
        bitstr3 = _BitString("0001")
        bitstr4 = _BitString("0010")
        self.assertRaises(TypeError, _BitString, "10O1")
        self.assertEqual(bitstr1 & bitstr2, _BitString("0001"))
        self.assertEqual(bitstr1 | bitstr2, _BitString("0111"))
        self.assertEqual(bitstr1 ^ bitstr2, _BitString("0110"))
        self.assertFalse(bitstr1.contains(bitstr2))
        self.assertTrue(bitstr1.contains(bitstr1))
        self.assertTrue(bitstr1.contains(bitstr3))
        self.assertTrue(bitstr1.contains(bitstr4))
        self.assertFalse(bitstr1.independent(bitstr2))
        self.assertFalse(bitstr1.independent(bitstr4))
        self.assertTrue(bitstr2.independent(bitstr4))
        self.assertTrue(bitstr3.independent(bitstr4))
        self.assertFalse(bitstr1.iscompatible(bitstr2))
        self.assertTrue(bitstr1.iscompatible(bitstr3))
        self.assertTrue(bitstr1.iscompatible(bitstr4))
        self.assertTrue(bitstr2.iscompatible(bitstr4))
        self.assertTrue(bitstr3.iscompatible(bitstr4))


class ConsensusTest(unittest.TestCase):
    """Test for consensus methods."""

    def setUp(self):
        self.trees = list(
            Phylo.parse(support.DATA / "TreeConstruction" / "trees.tre", "newick")
        )

    def test_count_clades(self):
        bitstr_counts, len_trees = Consensus._count_clades(self.trees)
        self.assertEqual(len_trees, len(self.trees))
        self.assertEqual(len(bitstr_counts), 6)
        # Bits stand for Alpha, Beta, Delta, Epsilon and Gamma, in that order
        self.assertEqual(bitstr_counts[_BitString("11111")][0], 3)
        self.assertEqual(bitstr_counts[_BitString("00110")][0], 2)
        self.assertEqual(bitstr_counts[_BitString("11001")][0], 3)
        self.assertEqual(bitstr_counts[_BitString("01001")][0], 2)
        self.assertEqual(bitstr_counts[_BitString("11000")][0], 1)
        self.assertEqual(bitstr_counts[_BitString("11101")][0], 1)

    def test_strict_consensus(self):
        ref_trees = list(
            Phylo.parse(support.DATA / "TreeConstruction" / "strict_refs.tre", "newick")
        )
        # three trees
        consensus_tree = Consensus.strict_consensus(self.trees)
        # tree_file = StringIO()
        # Phylo.write(consensus_tree, tree_file, 'newick')
        self.assertTrue(Consensus._equal_topology(consensus_tree, ref_trees[0]))
        # tree 1 and tree 2
        consensus_tree = Consensus.strict_consensus(self.trees[:2])
        # tree_file = StringIO()
        # Phylo.write(consensus_tree, tree_file, 'newick')
        self.assertTrue(Consensus._equal_topology(consensus_tree, ref_trees[1]))
        # tree 1 and tree 3
        consensus_tree = Consensus.strict_consensus(self.trees[::2])
        # tree_file = StringIO()
        # Phylo.write(consensus_tree, tree_file, 'newick')
        self.assertTrue(Consensus._equal_topology(consensus_tree, ref_trees[2]))
        # tree_file.close()

    def test_majority_consensus(self):
        ref_trees = Phylo.parse(
            support.DATA / "TreeConstruction" / "majority_ref.tre", "newick"
        )
        ref_tree = next(ref_trees)
        consensus_tree = Consensus.majority_consensus(self.trees)
        self.assertTrue(Consensus._equal_topology(consensus_tree, ref_tree))
        ref_tree = next(ref_trees)
        consensus_tree = Consensus.majority_consensus(self.trees, 1)
        self.assertTrue(Consensus._equal_topology(consensus_tree, ref_tree))

    def test_majority_consensus_identical_trees(self):
        # (D,E) used to be left out: once (C,(A,B)) was placed, the root had
        # three children, and that ended the search.
        trees = [read_newick("(((A,B),C),(D,E));") for _ in range(3)]
        for cutoff in (0, 0.5, 1):
            with self.subTest(cutoff=cutoff):
                consensus_tree = Consensus.majority_consensus(trees, cutoff)
                self.assertEqual(
                    inner_clades(consensus_tree),
                    {
                        frozenset("AB"): 100.0,
                        frozenset("ABC"): 100.0,
                        frozenset("DE"): 100.0,
                    },
                )

    def test_majority_consensus_three_way_roots(self):
        # Trees drawn unrooted, with three clades at the root, as nj draws
        # them. Once (A,B) and (D,E) are placed, the consensus is fully
        # resolved as an unrooted tree. (A,B,C) is the same split as (D,E),
        # so it is not added, though one tree draws its tips as a clade.
        trees = [
            read_newick("((A,B),C,(D,E));"),
            read_newick("((A,B),C,(D,E));"),
            read_newick("(((A,B),C),D,E);"),
        ]
        consensus_tree = Consensus.majority_consensus(trees)
        self.assertEqual(len(consensus_tree.root.clades), 3)
        self.assertEqual(
            inner_clades(consensus_tree).keys(), {frozenset("AB"), frozenset("DE")}
        )

    def test_adam_consensus(self):
        ref_trees = list(
            Phylo.parse(support.DATA / "TreeConstruction" / "adam_refs.tre", "newick")
        )
        # three trees
        consensus_tree = Consensus.adam_consensus(self.trees)
        # tree_file = '/home/yeyanbo/adam.tres'
        # tree_file = StringIO()
        # Phylo.write(consensus_tree, tree_file, 'newick')
        self.assertTrue(Consensus._equal_topology(consensus_tree, ref_trees[0]))
        # tree 1 and tree 2
        consensus_tree = Consensus.adam_consensus(self.trees[:2])
        # tree_file = StringIO()
        # Phylo.write(consensus_tree, tree_file, 'newick')
        self.assertTrue(Consensus._equal_topology(consensus_tree, ref_trees[1]))
        # tree 1 and tree 3
        consensus_tree = Consensus.adam_consensus(self.trees[::2])
        # tree_file = StringIO()
        # Phylo.write(consensus_tree, tree_file, 'newick')
        self.assertTrue(Consensus._equal_topology(consensus_tree, ref_trees[2]))
        # tree_file.close()

    def test_get_support(self):
        support_tree = Consensus.get_support(self.trees[0], self.trees)
        clade = support_tree.common_ancestor(
            [support_tree.find_any(name="Beta"), support_tree.find_any(name="Gamma")]
        )
        self.assertEqual(clade.confidence, 2 * 100.0 / 3)
        clade = support_tree.common_ancestor(
            [support_tree.find_any(name="Alpha"), support_tree.find_any(name="Beta")]
        )
        self.assertEqual(clade.confidence, 3 * 100.0 / 3)
        clade = support_tree.common_ancestor(
            [support_tree.find_any(name="Delta"), support_tree.find_any(name="Epsilon")]
        )
        self.assertEqual(clade.confidence, 2 * 100.0 / 3)


class TipOrderTest(unittest.TestCase):
    """Consensus methods must not depend on the order a tree lists its tips in."""

    def setUp(self):
        # One topology, ((A,B),C) and (D,E), listing its tips in three orders
        self.trees = [
            read_newick("(((A,B),C),(D,E));"),
            read_newick("((D,E),(C,(A,B)));"),
            read_newick("((E,D),((B,A),C));"),
        ]
        self.clades = {
            frozenset("AB"): 100.0,
            frozenset("ABC"): 100.0,
            frozenset("DE"): 100.0,
        }

    def test_count_clades(self):
        bitstr_counts, tree_count = Consensus._count_clades(self.trees)
        self.assertEqual(tree_count, 3)
        # The root, AB, ABC and DE, each in all three trees
        counts = sorted(count for count, branch_length in bitstr_counts.values())
        self.assertEqual(counts, [3, 3, 3, 3])

    def test_strict_consensus(self):
        consensus_tree = Consensus.strict_consensus(self.trees)
        self.assertEqual(inner_clades(consensus_tree).keys(), self.clades.keys())

    def test_majority_consensus(self):
        for cutoff in (0, 0.5, 1):
            with self.subTest(cutoff=cutoff):
                consensus_tree = Consensus.majority_consensus(self.trees, cutoff)
                self.assertEqual(inner_clades(consensus_tree), self.clades)

    def test_majority_consensus_ties(self):
        # Each clade is in one tree of two. With cutoff 0 the first clade
        # placed wins, and ties go by tip names, not by the order in which
        # the first tree lists its tips.
        expected = {frozenset("AB"): 50.0, frozenset("CD"): 50.0}
        for newicks in (
            ["((A,B),(C,D));", "((A,C),(B,D));"],
            ["((A,C),(B,D));", "((A,B),(C,D));"],
            ["((D,B),(C,A));", "((B,A),(D,C));"],
        ):
            with self.subTest(newicks=newicks):
                trees = [read_newick(newick) for newick in newicks]
                consensus_tree = Consensus.majority_consensus(trees)
                self.assertEqual(inner_clades(consensus_tree), expected)

    def test_equal_topology(self):
        self.assertTrue(Consensus._equal_topology(self.trees[0], self.trees[1]))
        # Different trees: DE sits where AB does, in each tree's own tip order
        self.assertFalse(
            Consensus._equal_topology(
                read_newick("((A,B),C,D,E);"), read_newick("((D,E),A,B,C);")
            )
        )

    def test_unnamed_terminal(self):
        # Sorting the names must not fail on a tip without one (None)
        trees = [read_newick("(((A,B),C),(D,));"), read_newick("((,D),(C,(B,A)));")]
        consensus_tree = Consensus.majority_consensus(trees)
        self.assertEqual(
            inner_clades(consensus_tree),
            {
                frozenset("AB"): 100.0,
                frozenset("ABC"): 100.0,
                frozenset({"D", None}): 100.0,
            },
        )

    def test_unnamed_terminal_and_tip_named_none(self):
        # An unnamed tip (None) must not sort as a tie with a tip named "None",
        # or the two trees below, one topology in two tip orders, would list
        # their names in different orders and be taken for different taxa.
        trees = [read_newick("((A,None),(B,));"), read_newick("((,B),(None,A));")]
        expected = {frozenset({"A", "None"}), frozenset({"B", None})}
        consensus_tree = Consensus.strict_consensus(trees)
        self.assertEqual(inner_clades(consensus_tree).keys(), expected)
        consensus_tree = Consensus.majority_consensus(trees)
        self.assertEqual(inner_clades(consensus_tree), dict.fromkeys(expected, 100.0))

    def test_inconsistent_taxa(self):
        for newick in ("(((A,B),C),(D,F));", "(((A,B),C),D);", "((A,B),C,D,E,F);"):
            trees = [self.trees[0], read_newick(newick)]
            for consensus in (
                Consensus.strict_consensus,
                Consensus.majority_consensus,
            ):
                with self.subTest(newick=newick, consensus=consensus.__name__):
                    with self.assertRaisesRegex(ValueError, "Taxons"):
                        consensus(trees)


class BootstrapTest(unittest.TestCase):
    """Test for bootstrap methods."""

    def setUp(self):
        self.msa = AlignIO.read(support.DATA / "TreeConstruction" / "msa.phy", "phylip")
        self.alignment = Align.read(
            support.DATA / "TreeConstruction" / "msa.phy", "phylip"
        )

    def test_bootstrap(self):
        msa_list = list(Consensus.bootstrap(self.msa, 100))
        self.assertEqual(len(msa_list), 100)
        self.assertEqual(len(msa_list[0]), len(self.msa))
        self.assertEqual(len(msa_list[0][0]), len(self.msa[0]))

    def test_bootstrap_trees_msa(self):
        calculator = DistanceCalculator("blosum62")
        constructor = DistanceTreeConstructor(calculator)
        trees = list(Consensus.bootstrap_trees(self.msa, 100, constructor))
        self.assertEqual(len(trees), 100)
        self.assertIsInstance(trees[0], BaseTree.Tree)

    def test_bootstrap_trees(self):
        calculator = DistanceCalculator("blosum62")
        constructor = DistanceTreeConstructor(calculator)
        trees = list(Consensus.bootstrap_trees(self.alignment, 100, constructor))
        self.assertEqual(len(trees), 100)
        self.assertIsInstance(trees[0], BaseTree.Tree)

    def test_majority_consensus_support(self):
        calculator = DistanceCalculator("identity")
        constructor = DistanceTreeConstructor(calculator, "upgma")
        state = random.getstate()
        random.seed(1)
        try:
            trees = list(Consensus.bootstrap_trees(self.alignment, 30, constructor))
        finally:
            random.setstate(state)
        # The replicate trees list their tips in different orders, as trees
        # built from bootstrap replicates usually do; that is the case tested.
        tip_orders = {
            tuple(term.name for term in tree.get_terminals()) for tree in trees
        }
        self.assertGreater(len(tip_orders), 1)
        replicate_clades = [inner_clades(tree).keys() for tree in trees]
        for cutoff in (0, 0.5):
            with self.subTest(cutoff=cutoff):
                consensus_tree = Consensus.majority_consensus(trees, cutoff)
                clades = inner_clades(consensus_tree)
                # Each clade's support is the share of replicates it is in
                for names, confidence in clades.items():
                    count = sum(names in replicate for replicate in replicate_clades)
                    self.assertAlmostEqual(confidence, 100.0 * count / len(trees))
                # and no clade in more than half the replicates is missing
                for replicate in replicate_clades:
                    for names in replicate:
                        count = sum(names in other for other in replicate_clades)
                        if 2 * count > len(trees):
                            self.assertIn(names, clades)

    def test_majority_consensus_nj_splits(self):
        # nj draws its unrooted trees with three clades at the root. Their
        # consensus must not add a clade that only splits the tips as one
        # already there does (the complement of a clade or of a tip): that
        # would add a root, not a split, with support counting only the
        # replicates that drew those tips as a clade.
        calculator = DistanceCalculator("identity")
        constructor = DistanceTreeConstructor(calculator, "nj")
        state = random.getstate()
        random.seed(1)
        try:
            trees = list(Consensus.bootstrap_trees(self.alignment, 50, constructor))
        finally:
            random.setstate(state)
        self.assertEqual({len(tree.root.clades) for tree in trees}, {3})
        names = frozenset(term.name for term in trees[0].get_terminals())
        consensus_tree = Consensus.majority_consensus(trees)
        clades = inner_clades(consensus_tree).keys()
        for clade in clades:
            complement = names - clade
            self.assertGreater(len(complement), 1)
            self.assertNotIn(complement, clades)

    def test_bootstrap_consensus_msa(self):
        calculator = DistanceCalculator("blosum62")
        constructor = DistanceTreeConstructor(calculator, "nj")
        tree = Consensus.bootstrap_consensus(
            self.msa, 100, constructor, Consensus.majority_consensus
        )
        self.assertIsInstance(tree, BaseTree.Tree)
        Phylo.write(tree, os.path.join(temp_dir, "bootstrap_consensus.tre"), "newick")

    def test_bootstrap_consensus(self):
        calculator = DistanceCalculator("blosum62")
        constructor = DistanceTreeConstructor(calculator, "nj")
        tree = Consensus.bootstrap_consensus(
            self.alignment, 100, constructor, Consensus.majority_consensus
        )
        self.assertIsInstance(tree, BaseTree.Tree)
        Phylo.write(tree, os.path.join(temp_dir, "bootstrap_consensus.tre"), "newick")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
