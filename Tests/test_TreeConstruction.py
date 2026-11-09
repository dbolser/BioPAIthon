# Copyright (C) 2013 by Yanbo Ye (yeyanbo289@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the Bio.Phylo.TreeConstruction module."""

import os
import tempfile
import unittest
from io import StringIO
from itertools import combinations

import support

from Bio import Align
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo import BaseTree
from Bio.Phylo import Consensus
from Bio.Phylo import TreeConstruction
from Bio.Phylo.TreeConstruction import _Matrix
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import NNITreeSearcher
from Bio.Phylo.TreeConstruction import ParsimonyScorer
from Bio.Phylo.TreeConstruction import ParsimonyTreeConstructor

temp_dir = tempfile.mkdtemp()


class DistanceMatrixTest(unittest.TestCase):
    """Test for DistanceMatrix construction and manipulation."""

    def setUp(self):
        self.names = ["Alpha", "Beta", "Gamma", "Delta"]
        self.matrix = [[0], [1, 0], [2, 3, 0], [4, 5, 6, 0]]

    def test_good_construction(self):
        dm = DistanceMatrix(self.names, self.matrix)
        self.assertIsInstance(dm, TreeConstruction.DistanceMatrix)
        self.assertEqual(dm.names[0], "Alpha")
        self.assertEqual(dm.matrix[2][1], 3)
        self.assertEqual(len(dm), 4)
        self.assertEqual(
            repr(dm),
            "DistanceMatrix(names=['Alpha', 'Beta', 'Gamma', 'Delta'], "
            "matrix=[[0], [1, 0], [2, 3, 0], [4, 5, 6, 0]])",
        )

    def test_bad_construction(self):
        self.assertRaises(
            TypeError,
            DistanceMatrix,
            ["Alpha", 100, "Gamma", "Delta"],
            [[0], [0.1, 0], [0.2, 0.3, 0], [0.4, 0.5, 0.6, 0]],
        )
        self.assertRaises(
            TypeError,
            DistanceMatrix,
            ["Alpha", "Beta", "Gamma", "Delta"],
            [[0], ["a"], [0.2, 0.3], [0.4, 0.5, 0.6]],
        )
        self.assertRaises(
            ValueError,
            DistanceMatrix,
            ["Alpha", "Alpha", "Gamma", "Delta"],
            [[0], [0.1], [0.2, 0.3], [0.4, 0.5, 0.6]],
        )
        self.assertRaises(
            ValueError,
            DistanceMatrix,
            ["Alpha", "Beta", "Gamma", "Delta"],
            [[0], [0.2, 0], [0.4, 0.5, 0.6]],
        )
        self.assertRaises(
            ValueError,
            DistanceMatrix,
            ["Alpha", "Beta", "Gamma", "Delta"],
            [[0], [0.1], [0.2, 0.3, 0.4], [0.4, 0.5, 0.6]],
        )

    def test_good_manipulation(self):
        dm = DistanceMatrix(self.names, self.matrix)
        # getitem
        self.assertEqual(dm[1], [1, 0, 3, 5])
        self.assertEqual(dm[2, 1], 3)
        self.assertEqual(dm[2][1], 3)
        self.assertEqual(dm[1, 2], 3)
        self.assertEqual(dm[1][2], 3)
        self.assertEqual(dm["Alpha"], [0, 1, 2, 4])
        self.assertEqual(dm["Gamma", "Delta"], 6)
        # setitem
        dm["Alpha"] = [0, 10, 20, 40]
        self.assertEqual(dm["Alpha"], [0, 10, 20, 40])
        # delitem insert item
        del dm[1]
        self.assertEqual(dm.names, ["Alpha", "Gamma", "Delta"])
        self.assertEqual(dm.matrix, [[0], [20, 0], [40, 6, 0]])
        dm.insert("Beta", [1, 0, 3, 5], 1)
        self.assertEqual(dm.names, self.names)
        self.assertEqual(dm.matrix, [[0], [1, 0], [20, 3, 0], [40, 5, 6, 0]])
        del dm["Alpha"]
        self.assertEqual(dm.names, ["Beta", "Gamma", "Delta"])
        self.assertEqual(dm.matrix, [[0], [3, 0], [5, 6, 0]])
        dm.insert("Alpha", [1, 2, 4, 0])
        self.assertEqual(dm.names, ["Beta", "Gamma", "Delta", "Alpha"])
        self.assertEqual(dm.matrix, [[0], [3, 0], [5, 6, 0], [1, 2, 4, 0]])

    def test_bad_manipulation(self):
        dm = DistanceMatrix(self.names, self.matrix)
        # getitem
        self.assertRaises(ValueError, dm.__getitem__, "A")
        self.assertRaises(ValueError, dm.__getitem__, ("Alpha", "A"))
        self.assertRaises(TypeError, dm.__getitem__, (1, "A"))
        self.assertRaises(TypeError, dm.__getitem__, (1, 1.2))
        self.assertRaises(IndexError, dm.__getitem__, 6)
        self.assertRaises(IndexError, dm.__getitem__, (10, 10))
        # setitem: item or index test
        self.assertRaises(ValueError, dm.__setitem__, "A", [1, 3, 4])
        self.assertRaises(ValueError, dm.__setitem__, ("Alpha", "A"), 4)
        self.assertRaises(TypeError, dm.__setitem__, (1, "A"), 3)
        self.assertRaises(TypeError, dm.__setitem__, (1, 1.2), 2)
        self.assertRaises(IndexError, dm.__setitem__, 6, [1, 3, 4])
        self.assertRaises(IndexError, dm.__setitem__, (10, 10), 1)
        # setitem: value test
        self.assertRaises(ValueError, dm.__setitem__, 0, [1, 2])
        self.assertRaises(TypeError, dm.__setitem__, ("Alpha", "Beta"), "a")
        self.assertRaises(TypeError, dm.__setitem__, "Alpha", ["a", "b", "c"])

    def test_format_phylip(self):
        dm = DistanceMatrix(self.names, self.matrix)
        handle = StringIO()
        dm.format_phylip(handle)
        lines = handle.getvalue().splitlines()
        self.assertEqual(len(lines), len(dm) + 1)
        self.assertTrue(lines[0].endswith(str(len(dm))))
        for name, line in zip(self.names, lines[1:]):
            self.assertTrue(line.startswith(name))


class DistanceCalculatorTest(unittest.TestCase):
    """Test DistanceCalculator."""

    def test_known_matrices_msa(self):
        msa = AlignIO.read(support.DATA / "TreeConstruction" / "msa.phy", "phylip")

        calculator = DistanceCalculator("identity")
        dm = calculator.get_distance(msa)
        self.assertEqual(dm["Alpha", "Beta"], 1 - 10 / 13)

        calculator = DistanceCalculator("blastn")
        dm = calculator.get_distance(msa)
        self.assertEqual(dm["Alpha", "Beta"], 1 - 38 / 65)

        calculator = DistanceCalculator("trans")
        dm = calculator.get_distance(msa)
        self.assertEqual(dm["Alpha", "Beta"], 1 - 54 / 65)

        calculator = DistanceCalculator("blosum62")
        dm = calculator.get_distance(msa)
        self.assertEqual(dm["Alpha", "Beta"], 1 - 53 / 84)

    def test_known_matrices(self):
        aln = Align.read(support.DATA / "TreeConstruction" / "msa.phy", "phylip")

        calculator = DistanceCalculator("identity")
        dm = calculator.get_distance(aln)
        self.assertEqual(dm["Alpha", "Beta"], 1 - 10 / 13)

        calculator = DistanceCalculator("blastn")
        dm = calculator.get_distance(aln)
        self.assertEqual(dm["Alpha", "Beta"], 1 - 38 / 65)

        calculator = DistanceCalculator("trans")
        dm = calculator.get_distance(aln)
        self.assertEqual(dm["Alpha", "Beta"], 1 - 54 / 65)

        calculator = DistanceCalculator("blosum62")
        dm = calculator.get_distance(aln)
        self.assertEqual(dm["Alpha", "Beta"], 1 - 53 / 84)

    def test_nonmatching_seqs_msa(self):
        aln = AlignIO.read(StringIO(">Alpha\nA-A--\n>Gamma\n-Y-Y-"), "fasta")
        # With a proper scoring matrix -- no matches
        dmat = DistanceCalculator("blosum62").get_distance(aln)
        self.assertEqual(dmat["Alpha", "Alpha"], 0.0)
        self.assertEqual(dmat["Alpha", "Gamma"], 1.0)
        # Comparing characters only -- 4 misses, 1 match
        dmat = DistanceCalculator().get_distance(aln)
        self.assertEqual(dmat["Alpha", "Alpha"], 0.0)
        self.assertAlmostEqual(dmat["Alpha", "Gamma"], 4.0 / 5.0)

    def test_nonmatching_seqs(self):
        aln = Align.read(
            StringIO(">Alpha\nA-A--\n>Beta\nXXXXX\n>Gamma\n-Y-Y-"), "fasta"
        )
        # With a proper scoring matrix -- no matches
        dmat = DistanceCalculator("blosum62").get_distance(aln)
        self.assertEqual(dmat["Alpha", "Alpha"], 0.0)
        self.assertEqual(dmat["Alpha", "Gamma"], 1.0)
        # Comparing characters only -- 4 misses, 1 match
        dmat = DistanceCalculator().get_distance(aln)
        self.assertEqual(dmat["Alpha", "Alpha"], 0.0)
        self.assertAlmostEqual(dmat["Alpha", "Gamma"], 4.0 / 5.0)

    def test_identity_skip_letters(self):
        # Identical sequences containing gap columns; skipped positions
        # must not inflate the distance (biopython/biopython#5332).
        aln = AlignIO.read(StringIO(">Alpha\nAC--GT\n>Gamma\nAC--GT"), "fasta")
        calculator = DistanceCalculator("identity", skip_letters=("-", "*"))
        self.assertEqual(calculator.get_distance(aln)["Alpha", "Gamma"], 0.0)
        # A gap column should count the same as no gap column.
        ungapped = AlignIO.read(StringIO(">Alpha\nACGT\n>Gamma\nACGA"), "fasta")
        gapped = AlignIO.read(StringIO(">Alpha\nAC--GT\n>Gamma\nAC--GA"), "fasta")
        self.assertEqual(
            calculator.get_distance(ungapped)["Alpha", "Gamma"],
            calculator.get_distance(gapped)["Alpha", "Gamma"],
        )

    def test_identity_skip_letters_gap_in_one_row(self):
        # A column with a gap in only one of the two rows is skipped for that
        # pair, with an Alignment as well as a MultipleSeqAlignment. Align.read
        # drops the columns that are gaps in every row, so the identical rows
        # above would not exercise skipping with an Alignment.
        calculator = DistanceCalculator("identity", skip_letters=("-", "*"))
        for fasta, distance in (
            (">Alpha\nAC--GT\n>Gamma\nACTTGT\n", 0.0),  # was 1 - 4/6
            (">Alpha\nAC--GA\n>Gamma\nACTTGT\n", 0.25),  # was 1 - 3/6
        ):
            for aln in (
                AlignIO.read(StringIO(fasta), "fasta"),
                Align.read(StringIO(fasta), "fasta"),
            ):
                with self.subTest(fasta=fasta, type=type(aln).__name__):
                    self.assertEqual(len(aln[0]), 6)
                    dm = calculator.get_distance(aln)
                    self.assertEqual(dm["Alpha", "Gamma"], distance)


class DistanceTreeConstructorTest(unittest.TestCase):
    """Test DistanceTreeConstructor."""

    def setUp(self):
        self.msa = AlignIO.read(support.DATA / "TreeConstruction" / "msa.phy", "phylip")
        calculator = DistanceCalculator("blosum62")
        self.dm_msa = calculator.get_distance(self.msa)
        self.constructor_msa = DistanceTreeConstructor(calculator)
        self.alignment = Align.read(
            support.DATA / "TreeConstruction" / "msa.phy", "phylip"
        )
        calculator = DistanceCalculator("blosum62")
        self.dm = calculator.get_distance(self.alignment)
        self.constructor = DistanceTreeConstructor(calculator)
        # UPGMA and WPGMA give different trees for this matrix
        self.dm_abcde = DistanceMatrix(
            ["A", "B", "C", "D", "E"],
            [[0], [2, 0], [4, 4, 0], [5, 5, 11, 0], [10, 10, 10, 7.5, 0]],
        )

    def test_upgma_msa(self):
        tree = self.constructor.upgma(self.dm_msa)
        self.assertIsInstance(tree, BaseTree.Tree)
        # tree_file = StringIO()
        # Phylo.write(tree, tree_file, 'newick')
        ref_tree = Phylo.read(support.DATA / "TreeConstruction" / "upgma.tre", "newick")
        self.assertTrue(Consensus._equal_topology(tree, ref_tree))
        # ref_tree.close()

    def test_upgma(self):
        tree = self.constructor.upgma(self.dm)
        self.assertIsInstance(tree, BaseTree.Tree)
        # tree_file = StringIO()
        # Phylo.write(tree, tree_file, 'newick')
        ref_tree = Phylo.read(support.DATA / "TreeConstruction" / "upgma.tre", "newick")
        self.assertTrue(Consensus._equal_topology(tree, ref_tree))
        # check for equal distance of all terminal nodes from the root
        ref_tree.root_at_midpoint()
        for len1, len2 in combinations(
            [depth for node, depth in ref_tree.depths().items() if node.is_terminal()],
            2,
        ):
            self.assertAlmostEqual(len1, len2)

    def check_tree_distances(self, tree, expected):
        """Check the path length between each pair of terminals, given as "AB"."""
        for pair, distance in expected.items():
            with self.subTest(pair=pair):
                self.assertAlmostEqual(tree.distance(*pair), distance)

    def test_upgma_weights_by_cluster_size(self):
        # biopython/biopython#1208: upgma() used to average the distances of
        # the two clusters it merged, which is WPGMA, so that D joined E at
        # 7.5 here instead of joining ABC at 7.
        tree = self.constructor.upgma(self.dm_abcde)
        # A and B join at 2, and C joins them at 4. The mean distance from D to
        # A, B and C is (5 + 5 + 11) / 3 = 7, less than D to E (7.5), so D joins
        # them next, and E joins all four at (10 + 10 + 10 + 7.5) / 4 = 9.375.
        self.check_tree_distances(
            tree,
            {
                "AB": 2,
                "AC": 4,
                "BC": 4,
                "AD": 7,
                "BD": 7,
                "CD": 7,
                "AE": 9.375,
                "BE": 9.375,
                "CE": 9.375,
                "DE": 9.375,
            },
        )

    def test_wpgma(self):
        tree = self.constructor.wpgma(self.dm_abcde)
        # WPGMA puts D at (5 + 11) / 2 = 8 from the AB-C cluster, whatever the
        # cluster sizes, so D joins E at 7.5 first, and the two clusters join
        # at (10 + 8) / 2 = 9.
        self.check_tree_distances(
            tree,
            {
                "AB": 2,
                "AC": 4,
                "BC": 4,
                "DE": 7.5,
                "AD": 9,
                "BD": 9,
                "CD": 9,
                "AE": 9,
                "BE": 9,
                "CE": 9,
            },
        )
        for dm in (self.dm_msa, self.dm):
            tree = self.constructor.wpgma(dm)
            ref_tree = Phylo.read(
                support.DATA / "TreeConstruction" / "wpgma.tre", "newick"
            )
            self.assertTrue(Consensus._equal_topology(tree, ref_tree))

    def test_build_tree_upgma_wpgma(self):
        calculator = DistanceCalculator("blosum62")
        for method in ("upgma", "wpgma"):
            with self.subTest(method=method):
                constructor = DistanceTreeConstructor(calculator, method)
                tree = constructor.build_tree(self.alignment)
                ref_tree = Phylo.read(
                    support.DATA / "TreeConstruction" / f"{method}.tre", "newick"
                )
                self.assertTrue(Consensus._equal_topology(tree, ref_tree))

    def test_nj_msa(self):
        tree = self.constructor.nj(self.dm_msa)
        self.assertIsInstance(tree, BaseTree.Tree)
        # tree_file = StringIO()
        # Phylo.write(tree, tree_file, 'newick')
        ref_tree = Phylo.read(support.DATA / "TreeConstruction" / "nj.tre", "newick")
        self.assertTrue(Consensus._equal_topology(tree, ref_tree))
        # ref_tree.close()

        # create a matrix of length 2
        calculator = DistanceCalculator("blosum62")
        self.min_dm = calculator.get_distance(self.msa)
        for i in range(len(self.min_dm) - 2):
            del self.min_dm[len(self.min_dm) - 1]

        min_tree = self.constructor.nj(self.min_dm)
        self.assertIsInstance(min_tree, BaseTree.Tree)

        ref_min_tree = Phylo.read(
            support.DATA / "TreeConstruction" / "nj_min.tre", "newick"
        )
        self.assertTrue(Consensus._equal_topology(min_tree, ref_min_tree))

    def test_nj(self):
        tree = self.constructor.nj(self.dm)
        self.assertIsInstance(tree, BaseTree.Tree)
        # tree_file = StringIO()
        # Phylo.write(tree, tree_file, 'newick')
        ref_tree = Phylo.read(support.DATA / "TreeConstruction" / "nj.tre", "newick")
        self.assertTrue(Consensus._equal_topology(tree, ref_tree))
        # ref_tree.close()

        # create a matrix of length 2
        calculator = DistanceCalculator("blosum62")
        self.min_dm = calculator.get_distance(self.msa)
        for i in range(len(self.min_dm) - 2):
            del self.min_dm[len(self.min_dm) - 1]

        min_tree = self.constructor.nj(self.min_dm)
        self.assertIsInstance(min_tree, BaseTree.Tree)

        ref_min_tree = Phylo.read(
            support.DATA / "TreeConstruction" / "nj_min.tre", "newick"
        )
        self.assertTrue(Consensus._equal_topology(min_tree, ref_min_tree))

    def test_built_tree_msa(self):
        tree = self.constructor.build_tree(self.msa)
        self.assertIsInstance(tree, BaseTree.Tree)
        # tree_file = StringIO()
        # Phylo.write(tree, tree_file, 'newick')
        ref_tree = Phylo.read(support.DATA / "TreeConstruction" / "nj.tre", "newick")
        self.assertTrue(Consensus._equal_topology(tree, ref_tree))
        # ref_tree.close()

    def test_built_tree(self):
        tree = self.constructor.build_tree(self.alignment)
        self.assertIsInstance(tree, BaseTree.Tree)
        # tree_file = StringIO()
        # Phylo.write(tree, tree_file, 'newick')
        ref_tree = Phylo.read(support.DATA / "TreeConstruction" / "nj.tre", "newick")
        self.assertTrue(Consensus._equal_topology(tree, ref_tree))
        # ref_tree.close()


class ParsimonyScorerTest(unittest.TestCase):
    """Test ParsimonyScorer."""

    def test_get_score_msa(self):
        aln = AlignIO.read(support.DATA / "TreeConstruction" / "msa.phy", "phylip")
        tree = Phylo.read(support.DATA / "TreeConstruction" / "upgma.tre", "newick")
        scorer = ParsimonyScorer()
        score = scorer.get_score(tree, aln)
        self.assertEqual(score, 2 + 1 + 2 + 2 + 1 + 1 + 1 + 3)

        alphabet = ["A", "T", "C", "G"]
        step_matrix = [[0], [2.5, 0], [2.5, 1, 0], [1, 2.5, 2.5, 0]]
        matrix = _Matrix(alphabet, step_matrix)
        scorer = ParsimonyScorer(matrix)
        score = scorer.get_score(tree, aln)
        self.assertEqual(score, 3.5 + 2.5 + 3.5 + 3.5 + 2.5 + 1 + 2.5 + 4.5)

        alphabet = [
            "A",
            "C",
            "D",
            "E",
            "F",
            "G",
            "H",
            "I",
            "K",
            "L",
            "M",
            "N",
            "P",
            "Q",
            "R",
            "1",
            "2",
            "T",
            "V",
            "W",
            "Y",
            "*",
            "-",
        ]
        step_matrix = [
            [0],
            [2, 0],
            [1, 2, 0],
            [1, 2, 1, 0],
            [2, 1, 2, 2, 0],
            [1, 1, 1, 1, 2, 0],
            [2, 2, 1, 2, 2, 2, 0],
            [2, 2, 2, 2, 1, 2, 2, 0],
            [2, 2, 2, 1, 2, 2, 2, 1, 0],
            [2, 2, 2, 2, 1, 2, 1, 1, 2, 0],
            [2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 0],
            [2, 2, 1, 2, 2, 2, 1, 1, 1, 2, 2, 0],
            [1, 2, 2, 2, 2, 2, 1, 2, 2, 1, 2, 2, 0],
            [2, 2, 2, 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 0],
            [2, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 0],
            [1, 1, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 2, 2, 0],
            [2, 1, 2, 2, 2, 1, 2, 1, 2, 2, 2, 1, 2, 2, 1, 2, 0],
            [1, 2, 2, 2, 2, 2, 2, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 0],
            [1, 2, 1, 1, 1, 1, 2, 1, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2, 0],
            [2, 1, 2, 2, 2, 1, 2, 2, 2, 1, 2, 3, 2, 2, 1, 1, 2, 2, 2, 0],
            [2, 1, 1, 2, 1, 2, 1, 2, 2, 2, 3, 1, 2, 2, 2, 1, 2, 2, 2, 2, 0],
            [2, 1, 2, 1, 2, 1, 2, 2, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 0],
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0],
        ]

        matrix = _Matrix(alphabet, step_matrix)
        scorer = ParsimonyScorer(matrix)
        score = scorer.get_score(tree, aln)
        self.assertEqual(score, 3 + 1 + 3 + 3 + 2 + 1 + 2 + 5)

    def test_get_score(self):
        aln = Align.read(support.DATA / "TreeConstruction" / "msa.phy", "phylip")
        tree = Phylo.read(support.DATA / "TreeConstruction" / "upgma.tre", "newick")
        scorer = ParsimonyScorer()
        score = scorer.get_score(tree, aln)
        self.assertEqual(score, 2 + 1 + 2 + 2 + 1 + 1 + 1 + 3)

        alphabet = ["A", "T", "C", "G"]
        step_matrix = [[0], [2.5, 0], [2.5, 1, 0], [1, 2.5, 2.5, 0]]
        matrix = _Matrix(alphabet, step_matrix)
        scorer = ParsimonyScorer(matrix)
        score = scorer.get_score(tree, aln)
        self.assertEqual(score, 3.5 + 2.5 + 3.5 + 3.5 + 2.5 + 1 + 2.5 + 4.5)

        alphabet = [
            "A",
            "C",
            "D",
            "E",
            "F",
            "G",
            "H",
            "I",
            "K",
            "L",
            "M",
            "N",
            "P",
            "Q",
            "R",
            "1",
            "2",
            "T",
            "V",
            "W",
            "Y",
            "*",
            "-",
        ]
        step_matrix = [
            [0],
            [2, 0],
            [1, 2, 0],
            [1, 2, 1, 0],
            [2, 1, 2, 2, 0],
            [1, 1, 1, 1, 2, 0],
            [2, 2, 1, 2, 2, 2, 0],
            [2, 2, 2, 2, 1, 2, 2, 0],
            [2, 2, 2, 1, 2, 2, 2, 1, 0],
            [2, 2, 2, 2, 1, 2, 1, 1, 2, 0],
            [2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 0],
            [2, 2, 1, 2, 2, 2, 1, 1, 1, 2, 2, 0],
            [1, 2, 2, 2, 2, 2, 1, 2, 2, 1, 2, 2, 0],
            [2, 2, 2, 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 0],
            [2, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 0],
            [1, 1, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 2, 2, 0],
            [2, 1, 2, 2, 2, 1, 2, 1, 2, 2, 2, 1, 2, 2, 1, 2, 0],
            [1, 2, 2, 2, 2, 2, 2, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 0],
            [1, 2, 1, 1, 1, 1, 2, 1, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2, 0],
            [2, 1, 2, 2, 2, 1, 2, 2, 2, 1, 2, 3, 2, 2, 1, 1, 2, 2, 2, 0],
            [2, 1, 1, 2, 1, 2, 1, 2, 2, 2, 3, 1, 2, 2, 2, 1, 2, 2, 2, 2, 0],
            [2, 1, 2, 1, 2, 1, 2, 2, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 0],
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0],
        ]

        matrix = _Matrix(alphabet, step_matrix)
        scorer = ParsimonyScorer(matrix)
        score = scorer.get_score(tree, aln)
        self.assertEqual(score, 3 + 1 + 3 + 3 + 2 + 1 + 2 + 5)


class NNITreeSearcherTest(unittest.TestCase):
    """Test NNITreeSearcher."""

    def test_get_neighbors(self):
        tree = Phylo.read(support.DATA / "TreeConstruction" / "upgma.tre", "newick")
        alphabet = ["A", "T", "C", "G"]
        step_matrix = [[0], [2.5, 0], [2.5, 1, 0], [1, 2.5, 2.5, 0]]
        matrix = _Matrix(alphabet, step_matrix)
        scorer = ParsimonyScorer(matrix)
        searcher = NNITreeSearcher(scorer)
        trees = searcher._get_neighbors(tree)
        self.assertEqual(len(trees), 2 * (5 - 3))
        Phylo.write(trees, os.path.join(temp_dir, "neighbor_trees.tre"), "newick")


class ParsimonyTreeConstructorTest(unittest.TestCase):
    """Test ParsimonyTreeConstructor."""

    def test_build_tree_msa(self):
        aln = AlignIO.read(support.DATA / "TreeConstruction" / "msa.phy", "phylip")
        tree1 = Phylo.read(support.DATA / "TreeConstruction" / "upgma.tre", "newick")
        tree2 = Phylo.read(support.DATA / "TreeConstruction" / "nj.tre", "newick")
        alphabet = ["A", "T", "C", "G"]
        step_matrix = [[0], [2.5, 0], [2.5, 1, 0], [1, 2.5, 2.5, 0]]
        matrix = _Matrix(alphabet, step_matrix)
        scorer = ParsimonyScorer(matrix)
        searcher = NNITreeSearcher(scorer)
        constructor = ParsimonyTreeConstructor(searcher, tree1)
        best_tree = constructor.build_tree(aln)
        Phylo.write(best_tree, os.path.join(temp_dir, "pars1.tre"), "newick")
        constructor.starting_tree = tree2
        best_tree = constructor.build_tree(aln)
        Phylo.write(best_tree, os.path.join(temp_dir, "pars2.tre"), "newick")
        constructor.starting_tree = None
        best_tree = constructor.build_tree(aln)
        Phylo.write(best_tree, os.path.join(temp_dir, "pars3.tre"), "newick")

    def test_build_tree(self):
        aln = Align.read(support.DATA / "TreeConstruction" / "msa.phy", "phylip")
        tree1 = Phylo.read(support.DATA / "TreeConstruction" / "upgma.tre", "newick")
        tree2 = Phylo.read(support.DATA / "TreeConstruction" / "nj.tre", "newick")
        alphabet = ["A", "T", "C", "G"]
        step_matrix = [[0], [2.5, 0], [2.5, 1, 0], [1, 2.5, 2.5, 0]]
        matrix = _Matrix(alphabet, step_matrix)
        scorer = ParsimonyScorer(matrix)
        searcher = NNITreeSearcher(scorer)
        constructor = ParsimonyTreeConstructor(searcher, tree1)
        best_tree = constructor.build_tree(aln)
        Phylo.write(best_tree, os.path.join(temp_dir, "pars1.tre"), "newick")
        constructor.starting_tree = tree2
        best_tree = constructor.build_tree(aln)
        Phylo.write(best_tree, os.path.join(temp_dir, "pars2.tre"), "newick")
        constructor.starting_tree = None
        best_tree = constructor.build_tree(aln)
        Phylo.write(best_tree, os.path.join(temp_dir, "pars3.tre"), "newick")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
