# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the Bio.Phylo module."""

import math
import os
import random
import tempfile
import unittest
import warnings
from io import StringIO

import support

from Bio import BiopythonWarning
from Bio import Phylo
from Bio.Nexus.Nexus import NexusError
from Bio.Phylo import Newick
from Bio.Phylo import PhyloXML
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

# Example Newick and Nexus files
EX_NEWICK = support.DATA / "Nexus" / "int_node_labels.nwk"
EX_NEWICK2 = support.DATA / "Nexus" / "test.new"
EX_NEXUS = support.DATA / "Nexus" / "test_Nexus_input.nex"
EX_NEXUS2 = support.DATA / "Nexus" / "bats.nex"
EX_NEWICK_BOM = support.DATA / "Nexus" / "ByteOrderMarkFile.nwk"

# Example PhyloXML files
EX_APAF = support.DATA / "PhyloXML" / "apaf.xml"
EX_BCL2 = support.DATA / "PhyloXML" / "bcl_2.xml"
EX_DIST = support.DATA / "PhyloXML" / "distribution.xml"
EX_PHYLO = support.DATA / "PhyloXML" / "phyloxml_examples.xml"


class IOTests(unittest.TestCase):
    """Tests for parsing and writing the supported formats."""

    def test_newick_read_single1(self):
        """Read first Newick file with one tree."""
        tree = Phylo.read(EX_NEWICK, "newick")
        self.assertEqual(len(tree.get_terminals()), 28)

    def test_newick_read_single2(self):
        """Read second Newick file with one tree."""
        tree = Phylo.read(EX_NEWICK2, "newick")
        self.assertEqual(len(tree.get_terminals()), 33)
        self.assertEqual(tree.find_any("Homo sapiens").comment, "modern human")
        self.assertEqual(
            tree.find_any("Equus caballus").comment,
            "wild horse; also 'Equus ferus caballus'",
        )
        self.assertEqual(tree.root.confidence, 80)
        tree = Phylo.read(EX_NEWICK2, "newick", comments_are_confidence=True)
        self.assertEqual(tree.root.confidence, 100)

    def test_newick_read_single3(self):
        """Read Nexus file with one tree."""
        tree = Phylo.read(EX_NEXUS2, "nexus")
        self.assertEqual(len(tree.get_terminals()), 658)

    def test_unicode_exception(self):
        """Read a Newick file with a unicode byte order mark (BOM)."""
        with open(EX_NEWICK_BOM, encoding="utf-8") as handle:
            tree = Phylo.read(handle, "newick")
        self.assertEqual(len(tree.get_terminals()), 3)

    def test_newick_read_multiple(self):
        """Parse a Nexus file with multiple trees."""
        trees = list(Phylo.parse(EX_NEXUS, "nexus"))
        self.assertEqual(len(trees), 3)
        for tree in trees:
            self.assertEqual(len(tree.get_terminals()), 9)

    def test_nexus_read_deep_tree(self):
        """Read a Nexus tree nested deeper than the recursion limit."""
        depth = 5000
        newick = "t1"
        for i in range(2, depth + 1):
            newick = f"({newick},t{i})"
        handle = StringIO(f"#NEXUS\nbegin trees;\ntree deep = {newick};\nend;\n")
        tree = Phylo.read(handle, "nexus")
        # Walk down the spine with a loop: BaseTree's traversals recurse too
        clade = tree.root
        for level in range(depth - 1):
            self.assertEqual(len(clade.clades), 2)
            self.assertEqual(clade.clades[1].name, f"t{depth - level}")
            clade = clade.clades[0]
        self.assertEqual(clade.clades, [])
        self.assertEqual(clade.name, "t1")

    def test_newick_write_deep_tree(self):
        """Write a Newick tree nested deeper than the recursion limit."""
        depth = 5000
        newick = "t1:1"
        for i in range(2, depth + 1):
            newick = f"({newick},t{i}:1):1"
        tree = Phylo.read(StringIO(newick + ";"), "newick")
        handle = StringIO()
        Phylo.write(tree, handle, "newick")
        self.assertEqual(handle.getvalue(), newick + ";\n")
        handle = StringIO()
        Phylo.write(tree, handle, "newick", plain=True)
        self.assertEqual(handle.getvalue(), newick.replace(":1", "") + ";\n")

    def test_nexus_read_values(self):
        """Nexus trees keep the values Bio.Nexus.Trees gave them."""
        handle = StringIO(
            "#NEXUS\nbegin trees;\n"
            "translate 1 'Homo sapiens', 2 Pan, 3 Gorilla;\n"
            "tree t = [&R] [&W 0.5] ((1:0.1[&rate=1.5],2)95:0.3,3)80;\n"
            "end;\n"
        )
        tree = Phylo.read(handle, "nexus")
        self.assertEqual(tree.name, "t")
        self.assertTrue(tree.rooted)
        self.assertEqual(tree.weight, 0.5)
        inner, gorilla = tree.root.clades
        human, pan = inner.clades
        # Names from the translate table are NEXUS-quoted
        self.assertEqual(
            [c.name for c in (human, pan, gorilla)],
            ["'Homo sapiens'", "Pan", "Gorilla"],
        )
        # Confidences are floats
        self.assertIsInstance(inner.confidence, float)
        self.assertEqual(inner.confidence, 95.0)
        self.assertEqual(inner.branch_length, 0.3)
        # A bare number after a clade is its branch length
        self.assertEqual(tree.root.branch_length, 80.0)
        self.assertIsNone(tree.root.confidence)
        # A missing branch length is 0.0
        self.assertEqual(pan.branch_length, 0.0)
        self.assertEqual(gorilla.branch_length, 0.0)
        # Comments keep their brackets
        self.assertEqual(human.comment, "[&rate=1.5]")
        self.assertEqual(human.branch_length, 0.1)

    def test_nexus_read_negative_zero(self):
        """A bare -0.0 after a clade is a branch length of -0.0, as before."""
        handle = StringIO("#NEXUS\nbegin trees;\ntree t = ((a,b)-0.0,c);\nend;\n")
        inner = Phylo.read(handle, "nexus").root.clades[0]
        self.assertEqual(math.copysign(1, inner.branch_length), -1)

    def test_nexus_read_quoted_number(self):
        """A quoted internal label is a name, even if it looks like a number."""
        handle = StringIO(
            "#NEXUS\nbegin trees;\ntree t = ((a,b)'95',(c,d)'0.5':0.2);\nend;\n"
        )
        first, second = Phylo.read(handle, "nexus").root.clades
        self.assertEqual(first.name, "95")
        self.assertIsNone(first.confidence)
        self.assertEqual(first.branch_length, 0.0)
        self.assertEqual(second.name, "0.5")
        self.assertIsNone(second.confidence)
        self.assertEqual(second.branch_length, 0.2)

    def test_nexus_read_comment_after_colon(self):
        """A comment between the colon and branch length, as BEAST writes it."""
        handle = StringIO(
            "#NEXUS\nbegin trees;\ntranslate 1 A, 2 B, 3 C;\n"
            "tree t = ((1:[&rate=1.5]0.1,2)95:[&rate=0.5] 0.3,3);\nend;\n"
        )
        inner, c = Phylo.read(handle, "nexus").root.clades
        a, b = inner.clades
        self.assertEqual([a.name, b.name, c.name], ["A", "B", "C"])
        self.assertEqual(a.branch_length, 0.1)
        self.assertEqual(a.comment, "[&rate=1.5]")
        self.assertEqual(inner.confidence, 95.0)
        self.assertEqual(inner.branch_length, 0.3)
        self.assertEqual(inner.comment, "[&rate=0.5]")

    def test_nexus_read_comment_inside_number(self):
        """A comment inside a number does not split it, as NEXUS says."""
        handle = StringIO(
            "#NEXUS\nbegin trees;\n"
            "tree t = (a:-[&x]0.1,(b,c)9[&y]5:0.2,d:2[&z]e-1);\nend;\n"
        )
        a, inner, d = Phylo.read(handle, "nexus").root.clades
        self.assertEqual([a.name, d.name], ["a", "d"])
        self.assertEqual(a.branch_length, -0.1)
        self.assertEqual(a.comment, "[&x]")
        self.assertEqual(inner.confidence, 95.0)
        self.assertEqual(inner.branch_length, 0.2)
        self.assertEqual(inner.comment, "[&y]")
        self.assertEqual(d.branch_length, 0.2)
        self.assertEqual(d.comment, "[&z]")

    # A MrBayes 3.2 consensus tree: each clade has its posterior in a comment
    # before the branch length and branch length statistics in one after it
    MRBAYES_CONSENSUS = """\
#NEXUS
begin trees;
   translate
      1 A,
      2 B,
      3 C;
   tree con_50_majrule = [&U] \
(1[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob(percent)="100"]\
:1.0e-01[&length_mean=1.0e-01,length_95%HPD={9.0e-02,1.1e-01}],\
(2[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob(percent)="100"]\
:2.0e-01[&length_mean=2.0e-01,length_95%HPD={1.9e-01,2.1e-01}],\
3[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob(percent)="100"]\
:3.0e-01[&length_mean=3.0e-01,length_95%HPD={2.9e-01,3.1e-01}])\
[&prob=9.5e-01,prob_stddev=1.0e-02,prob(percent)="95"]\
:4.0e-01[&length_mean=4.0e-01,length_95%HPD={3.9e-01,4.1e-01}]);
end;
"""

    def test_nexus_read_two_comments(self):
        """A clade with two comments keeps both, in file order."""
        tree = Phylo.read(StringIO(self.MRBAYES_CONSENSUS), "nexus")
        a, inner = tree.root.clades
        b, c = inner.clades
        self.assertEqual([a.name, b.name, c.name], ["A", "B", "C"])
        self.assertEqual(
            [a.branch_length, inner.branch_length, b.branch_length, c.branch_length],
            [0.1, 0.4, 0.2, 0.3],
        )
        self.assertEqual(
            a.comment,
            '[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob(percent)="100"]'
            "[&length_mean=1.0e-01,length_95%HPD={9.0e-02,1.1e-01}]",
        )
        # The clade posterior, before the branch length, is not lost
        self.assertEqual(
            inner.comment,
            '[&prob=9.5e-01,prob_stddev=1.0e-02,prob(percent)="95"]'
            "[&length_mean=4.0e-01,length_95%HPD={3.9e-01,4.1e-01}]",
        )
        self.assertIsNone(tree.root.comment)

    def test_nexus_write_comments(self):
        """Comments read from a NEXUS tree are written so they read back."""
        tree = Phylo.read(StringIO(self.MRBAYES_CONSENSUS), "nexus")
        handle = StringIO()
        Phylo.write(tree, handle, "nexus")
        self.assertIn(
            "A:0.1"
            '[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob(percent)="100"]'
            "[&length_mean=1.0e-01,length_95%HPD={9.0e-02,1.1e-01}]",
            handle.getvalue(),
        )
        handle.seek(0)
        tree2 = Phylo.read(handle, "nexus")
        self.assertEqual(
            [(c.name, c.branch_length, c.comment) for c in tree2.find_clades()],
            [(c.name, c.branch_length, c.comment) for c in tree.find_clades()],
        )
        # One comment, as BEAST writes, reads back too
        tree = Phylo.read(
            StringIO("#NEXUS\nbegin trees;\ntree t = (a[&rate=1.5]:0.1,b);\nend;\n"),
            "nexus",
        )
        handle = StringIO()
        Phylo.write(tree, handle, "nexus")
        handle.seek(0)
        self.assertEqual(
            Phylo.read(handle, "nexus").root.clades[0].comment, "[&rate=1.5]"
        )
        # Newick has no such comments, so its writer still escapes brackets
        self.assertEqual(
            tree.format("newick").strip(), r"(a:0.1[\[&rate=1.5\]],b:0):0;"
        )

    def test_nexus_read_translate_error(self):
        """A failed TRANSLATE lookup names the first missing taxon in preorder."""
        handle = StringIO(
            "#NEXUS\nbegin trees;\ntranslate 1 a, 2 b;\n"
            "tree t = (1,(2,3),4);\nend;\n"
        )
        with self.assertRaisesRegex(NexusError, "Unable to substitute 3 "):
            Phylo.read(handle, "nexus")

    def test_newick_write(self):
        """Parse a Nexus file with multiple trees."""
        # Tree with internal node labels
        mem_file = StringIO()
        tree = Phylo.read(StringIO("(A,B,(C,D)E)F;"), "newick")
        Phylo.write(tree, mem_file, "newick")
        mem_file.seek(0)
        tree2 = Phylo.read(mem_file, "newick")
        # Sanity check
        self.assertEqual(tree2.count_terminals(), 4)
        # Check internal node labels were retained
        internal_names = {c.name for c in tree2.get_nonterminals() if c is not None}
        self.assertEqual(internal_names, {"E", "F"})

    def test_newick_write_unnamed_clade_confidence(self):
        """Write the confidence of an unnamed clade as its label."""
        tree = Phylo.read(StringIO("((C:0.16,D:0.14)98:0.46,A:0.11);"), "newick")
        newick = tree.format("newick").strip()
        self.assertEqual(newick, "((C:0.16,D:0.14)98.00:0.46,A:0.11):0;")
        inner = Phylo.read(StringIO(newick), "newick").root.clades[0]
        self.assertEqual([inner.name, inner.confidence], [None, 98.0])

    def test_newick_write_named_clade_confidence(self):
        """Write the confidence of a named clade in a comment, not its label."""
        # Written after the name, as for an unnamed clade, the confidence gave
        # (A:1,B:1)N10.90:0.5, which reads back as a clade named N10.90.
        tree = Phylo.read(StringIO("(A:1,B:1)N1:0.5;"), "newick")
        tree.root.confidence = 0.9
        newick = tree.format("newick").strip()
        self.assertEqual(newick, "(A:1,B:1)N1:0.5[0.90];")
        root = Phylo.read(StringIO(newick), "newick").root
        self.assertEqual(
            [root.name, root.confidence, root.comment], ["N1", None, "0.90"]
        )
        tree = Phylo.read(StringIO(newick), "newick", comments_are_confidence=True)
        self.assertEqual([tree.root.name, tree.root.confidence], ["N1", 0.9])
        # NEXUS ignores a plain comment, so only the name is read back
        handle = StringIO()
        Phylo.write(tree, handle, "nexus")
        handle.seek(0)
        root = Phylo.read(handle, "nexus").root
        self.assertEqual([root.name, root.confidence], ["N1", None])

    def test_newick_write_comments_are_confidence(self):
        """Write named clades read with comments_are_confidence=True as read."""
        newick = "(A:0.22,(B:0.12,(C:0.02,D:0.02)E:0.1[99])F:0.1[98]):0;"
        tree = Phylo.read(StringIO(newick), "newick", comments_are_confidence=True)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            self.assertEqual(tree.format("newick").strip(), newick)

    def test_newick_write_changed_confidence(self):
        """Write the new confidence of a named clade, not its old comment."""
        # Read with comments_are_confidence=True, a clade keeps its support
        # value in its comment too, which goes stale if the confidence changes
        newick = "(A:0.22,(B:0.12,(C:0.02,D:0.02)E:0.1[99])F:0.1[98]):0;"
        tree = Phylo.read(StringIO(newick), "newick", comments_are_confidence=True)
        for clade in tree.get_nonterminals():
            if clade.confidence is not None:
                clade.confidence /= 100
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            newick = tree.format("newick").strip()
        self.assertEqual(
            newick, "(A:0.22,(B:0.12,(C:0.02,D:0.02)E:0.1[0.99])F:0.1[0.98]):0;"
        )
        # Rerooting moves the support values of Q and R, but not their comments
        newick = "((((A:0.1,B:0.1)P:0.1[91],C:0.1)Q:0.1[82],D:0.1)R:0.1[73],G:0.1);"
        tree = Phylo.read(StringIO(newick), "newick", comments_are_confidence=True)
        tree.root_with_outgroup("A")
        newick = tree.format("newick").strip()
        tree = Phylo.read(StringIO(newick), "newick", comments_are_confidence=True)
        support = {c.name: c.confidence for c in tree.get_nonterminals() if c.name}
        self.assertEqual([support["Q"], support["R"]], [91, 82])

    def test_newick_write_named_clade_confidence_and_comment(self):
        """Write a named clade that has both a confidence and a comment."""
        root = Newick.Clade(
            branch_length=0.1,
            name="E",
            clades=[Newick.Clade(name="A"), Newick.Clade(name="B")],
            confidence=50,
            comment="note",
        )
        tree = Newick.Tree(root=root)
        with self.assertWarnsRegex(
            BiopythonWarning, "^Clade 'E' has a comment, so its confidence 50 was"
        ):
            newick = tree.format("newick").strip()
        self.assertEqual(newick, "(A:0,B:0)E:0.1[note];")
        # A comment that reads as the confidence holds it, so nothing is lost
        root.comment = "50"
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            newick = tree.format("newick").strip()
        self.assertEqual(newick, "(A:0,B:0)E:0.1[50];")
        # A comment that reads as another number may be data, so it is kept
        root.comment = "40"
        with self.assertWarnsRegex(BiopythonWarning, "^Clade 'E' has a comment"):
            newick = tree.format("newick").strip()
        self.assertEqual(newick, "(A:0,B:0)E:0.1[40];")

    def test_newick_write_named_clade_numeric_comment(self):
        """Keep a numeric comment that the confidence was not read from."""
        # The confidence is read from the value after the colon, not [2024]
        newick = "(A:1,B:1)N1:0.9[2024];"
        tree = Phylo.read(StringIO(newick), "newick", values_are_confidence=True)
        self.assertEqual([tree.root.confidence, tree.root.comment], [0.9, "2024"])
        with self.assertWarnsRegex(BiopythonWarning, "^Clade 'N1' has a comment"):
            newick = tree.format("newick").strip()
        self.assertEqual(newick, "(A:0,B:0)N1:0[2024];")
        # Read with the defaults, a comment is not a confidence, so one set
        # later does not replace it
        tree = Phylo.read(StringIO("(A:1,B:1)N1:0.9[12.5];"), "newick")
        tree.root.confidence = 0.75
        with self.assertWarnsRegex(BiopythonWarning, "^Clade 'N1' has a comment"):
            newick = tree.format("newick").strip()
        self.assertEqual(newick, "(A:1,B:1)N1:0.9[12.5];")

    def test_newick_read_scinot(self):
        """Parse Newick branch lengths in scientific notation."""
        tree = Phylo.read(StringIO("(foo:1e-1,bar:0.1)"), "newick")
        clade_a = tree.clade[0]
        self.assertEqual(clade_a.name, "foo")
        self.assertAlmostEqual(clade_a.branch_length, 0.1)

    def test_newick_read_whitespace_after_colon(self):
        """Read a branch length with any whitespace after the colon."""
        # Newick allows blanks, tabs and newlines anywhere except inside a
        # label or a branch length. Only one space after the colon used to
        # work: after two, or a tab, the length was read as the clade's name.
        tree = Phylo.read(
            StringIO("((a:  0.1,b:\t0.2)0.95 : \t0.3,c:\n    0.4);"), "newick"
        )
        inner, c = tree.root.clades
        a, b = inner.clades
        self.assertEqual([a.name, b.name, inner.name, c.name], ["a", "b", None, "c"])
        self.assertEqual(
            [a.branch_length, b.branch_length, inner.branch_length, c.branch_length],
            [0.1, 0.2, 0.3, 0.4],
        )
        self.assertEqual(inner.confidence, 0.95)
        # Any other whitespace \s matches, as it ends an unquoted label too.
        # float() does not skip "\x1c" to "\x1f", nor does str.strip() on PyPy 3.10.
        for space in "\x1c\x1f\x85\xa0\u3000":
            with self.subTest(space=hex(ord(space))):
                tree = Phylo.read(StringIO(f"(a:{space}0.1,b);"), "newick")
                a, b = tree.root.clades
                self.assertEqual([a.name, a.branch_length], ["a", 0.1])

    def test_newick_read_colon_inside_quotes_and_comments(self):
        """A colon and spaces inside a quoted label or comment are kept."""
        tree = Phylo.read(StringIO("('a:  1':0.1,b[x:  2]:0.2);"), "newick")
        a, b = tree.root.clades
        self.assertEqual([a.name, a.branch_length], ["a:  1", 0.1])
        self.assertEqual([b.name, b.comment, b.branch_length], ["b", "x:  2", 0.2])

    def test_phylo_read_extra(self):
        """Additional tests to check correct parsing."""
        tree = Phylo.read(StringIO("(A:1, B:-2, (C:3, D:4):-2)"), "newick")
        self.assertEqual(tree.distance("A"), 1)
        self.assertEqual(tree.distance("B"), -2)
        self.assertEqual(tree.distance("C"), 1)
        self.assertEqual(tree.distance("D"), 2)

        tree = Phylo.read(StringIO("((A:1, B:-2):-5, (C:3, D:4):-2)"), "newick")
        self.assertEqual(tree.distance("A"), -4)
        self.assertEqual(tree.distance("B"), -7)
        self.assertEqual(tree.distance("C"), 1)
        self.assertEqual(tree.distance("D"), 2)

        tree = Phylo.read(StringIO("((:1, B:-2):-5, (C:3, D:4):-2)"), "newick")
        distances = {-4.0: 1, -7.0: 1, 1: 1, 2: 1}
        for x in tree.get_terminals():
            entry = int(tree.distance(x))
            distances[entry] -= distances[entry]
            self.assertEqual(distances[entry], 0)

        tree = Phylo.read(StringIO("((:\n1\n,\n B:-2):-5, (C:3, D:4):-2);"), "newick")
        distances = {-4.0: 1, -7.0: 1, 1: 1, 2: 1}
        for x in tree.get_terminals():
            entry = int(tree.distance(x))
            distances[entry] -= distances[entry]
            self.assertEqual(distances[entry], 0)

    def test_format_branch_length(self):
        """Custom format string for Newick branch length serialization."""
        tree = Phylo.read(StringIO("A:0.1;"), "newick")
        mem_file = StringIO()
        Phylo.write(tree, mem_file, "newick", format_branch_length="%.0e")
        value = mem_file.getvalue().strip()
        self.assertTrue(value.startswith("A:"))
        self.assertTrue(value.endswith(";"))
        self.assertEqual(value[2:-1], "%.0e" % 0.1)

    def test_newick_writer_preserves_small_branch_lengths(self):
        """Small nonzero branch lengths should not be rounded to zero."""
        tree = Phylo.read(StringIO("(A:0.000001,B:0.000002);"), "newick")

        mem_file = StringIO()
        Phylo.write(tree, mem_file, "newick")

        newick = mem_file.getvalue()
        self.assertNotIn(":0.00000", newick)
        self.assertIn(":1e-06", newick)
        self.assertIn(":2e-06", newick)

    def test_convert(self):
        """Convert a tree between all supported formats."""
        mem_file_1 = StringIO()
        mem_file_2 = StringIO()
        mem_file_3 = StringIO()
        Phylo.convert(EX_NEWICK, "newick", mem_file_1, "nexus")
        mem_file_1.seek(0)
        Phylo.convert(mem_file_1, "nexus", mem_file_2, "phyloxml")
        mem_file_2.seek(0)
        Phylo.convert(mem_file_2, "phyloxml", mem_file_3, "newick")
        mem_file_3.seek(0)
        tree = Phylo.read(mem_file_3, "newick")
        self.assertEqual(len(tree.get_terminals()), 28)

    def test_convert_phyloxml_binary(self):
        """Try writing phyloxml to a binary handle; fail on Py3."""
        trees = Phylo.parse(
            support.DATA / "PhyloXML" / "phyloxml_examples.xml", "phyloxml"
        )
        with tempfile.NamedTemporaryFile(mode="wb") as out_handle:
            self.assertRaises(TypeError, Phylo.write, trees, out_handle, "phyloxml")

    def test_convert_phyloxml_text(self):
        """Write phyloxml to a text handle."""
        trees = Phylo.parse(
            support.DATA / "PhyloXML" / "phyloxml_examples.xml", "phyloxml"
        )
        with tempfile.NamedTemporaryFile(mode="w") as out_handle:
            count = Phylo.write(trees, out_handle, "phyloxml")
        self.assertEqual(14, count)

    def test_convert_phyloxml_filename(self):
        """Write phyloxml to a given filename."""
        trees = Phylo.parse(
            support.DATA / "PhyloXML" / "phyloxml_examples.xml", "phyloxml"
        )
        out_handle = tempfile.NamedTemporaryFile(mode="w", delete=False)
        out_handle.close()
        tmp_filename = out_handle.name
        try:
            count = Phylo.write(trees, tmp_filename, "phyloxml")
        finally:
            os.remove(tmp_filename)
        self.assertEqual(14, count)

    def test_convert_phyloxml_to_newick_branch_length_only(self):
        """Write phyloxml with bootstrap values to newick format using branch_length_only=True"""
        trees = Phylo.parse(EX_APAF, "phyloxml")
        out_handle = tempfile.NamedTemporaryFile(mode="w", delete=False)
        out_handle.close()
        tmp_filename = out_handle.name
        try:
            Phylo.write(trees, tmp_filename, "newick", branch_length_only=True)
        finally:
            os.remove(tmp_filename)

    def test_convert_phyloxml_to_newick_named_clade_confidence(self):
        """Keep the name and support of a PhyloXML clade written as Newick."""
        tree = list(Phylo.parse(EX_PHYLO, "phyloxml"))[2]
        self.assertEqual([tree.clade[0].name, tree.clade[0].confidence], ["AB", 89])
        handle = StringIO()
        Phylo.write(tree, handle, "newick")
        handle.seek(0)
        tree = Phylo.read(handle, "newick", comments_are_confidence=True)
        self.assertEqual([tree.clade[0].name, tree.clade[0].confidence], ["AB", 89])

    def test_int_labels(self):
        """Read newick formatted tree with numeric labels."""
        tree = Phylo.read(
            StringIO("(((0:0.1,1:0.1)0.99:0.1,2:0.1)0.98:0.0);"), "newick"
        )
        self.assertEqual({leaf.name for leaf in tree.get_terminals()}, {"0", "1", "2"})

    def test_newick_escaped_quote_parse(self):
        """Parse a newick tree with an escaped quote in a label."""
        newick_str = "('Node''Name':0.00000)Root:0.00000;"
        tree = Phylo.read(StringIO(newick_str), "newick")
        self.assertEqual(tree.get_terminals()[0].name, "Node'Name")

    def test_newick_escaped_quote_write(self):
        """Write a newick tree with an escaped quote in a label."""
        tree = Phylo.BaseTree.Tree(root=Phylo.BaseTree.Clade(name="Root"))
        clade_with_quote = Phylo.BaseTree.Clade(name="Node'Name")
        tree.root.clades.append(clade_with_quote)
        mem_file = StringIO()
        Phylo.write(tree, mem_file, "newick")
        mem_file.seek(0)
        self.assertEqual(mem_file.read(), "('Node''Name':0)Root:0;\n")


class TreeTests(unittest.TestCase):
    """Tests for methods on BaseTree.Tree objects."""

    def test_randomized(self):
        """Tree.randomized: generate a new randomized tree."""
        for N in (2, 5, 20):
            tree = Phylo.BaseTree.Tree.randomized(N)
            self.assertEqual(tree.count_terminals(), N)
            self.assertEqual(tree.total_branch_length(), (N - 1) * 2)
            tree = Phylo.BaseTree.Tree.randomized(N, branch_length=2.0)
            self.assertEqual(tree.total_branch_length(), (N - 1) * 4)
        tree = Phylo.BaseTree.Tree.randomized(5, branch_stdev=0.5)
        self.assertEqual(tree.count_terminals(), 5)

    def test_root_with_outgroup(self):
        """Tree.root_with_outgroup: reroot at a given clade."""
        # On a large realistic tree, at a deep internal node
        tree = Phylo.read(EX_APAF, "phyloxml")
        orig_num_tips = len(tree.get_terminals())
        orig_tree_len = tree.total_branch_length()
        tree.root_with_outgroup("19_NEMVE", "20_NEMVE")
        self.assertEqual(orig_num_tips, len(tree.get_terminals()))
        self.assertAlmostEqual(orig_tree_len, tree.total_branch_length())
        # Now, at an external node
        tree.root_with_outgroup("1_BRAFL")
        self.assertEqual(orig_num_tips, len(tree.get_terminals()))
        self.assertAlmostEqual(orig_tree_len, tree.total_branch_length())
        # Specifying outgroup branch length mustn't change the total tree size
        tree.root_with_outgroup("2_BRAFL", outgroup_branch_length=0.5)
        self.assertEqual(orig_num_tips, len(tree.get_terminals()))
        self.assertAlmostEqual(orig_tree_len, tree.total_branch_length())
        tree.root_with_outgroup("36_BRAFL", "37_BRAFL", outgroup_branch_length=0.5)
        self.assertEqual(orig_num_tips, len(tree.get_terminals()))
        self.assertAlmostEqual(orig_tree_len, tree.total_branch_length())
        # On small contrived trees, testing edge cases
        for small_nwk in (
            "(A,B,(C,D));",
            "((E,F),((G,H)),(I,J));",
            "((Q,R),(S,T),(U,V));",
            "(X,Y);",
        ):
            tree = Phylo.read(StringIO(small_nwk), "newick")
            orig_tree_len = tree.total_branch_length()
            for node in list(tree.find_clades()):
                tree.root_with_outgroup(node)
                self.assertAlmostEqual(orig_tree_len, tree.total_branch_length())
        # A non-zero outgroup branch length on a direct child of the root must
        # split the existing branch rather than be added to it.
        tree = Phylo.read(StringIO("(A:1.0,B:2.0);"), "newick")
        tree.root_with_outgroup("A", outgroup_branch_length=0.5)
        self.assertAlmostEqual(tree.total_branch_length(), 3.0)
        self.assertAlmostEqual(tree.distance("A", "B"), 3.0)
        tree = Phylo.read(StringIO("(A:1.0,(B:1.0,C:1.0):1.0);"), "newick")
        tree.root_with_outgroup("A", outgroup_branch_length=0.5)
        self.assertAlmostEqual(tree.total_branch_length(), 4.0)
        self.assertAlmostEqual(tree.distance("A", "B"), 3.0)

    def test_root_at_midpoint(self):
        """Tree.root_at_midpoint: reroot at the tree's midpoint."""
        for treefname, fmt in [
            (EX_APAF, "phyloxml"),
            (EX_BCL2, "phyloxml"),
            (EX_NEWICK, "newick"),
        ]:
            tree = Phylo.read(treefname, fmt)
            orig_tree_len = tree.total_branch_length()
            # Total branch length does not change
            tree.root_at_midpoint()
            self.assertAlmostEqual(orig_tree_len, tree.total_branch_length())
            # Root is bifurcating
            self.assertEqual(len(tree.root.clades), 2)
            # Deepest tips under each child of the root are equally deep
            deep_dist_0 = max(tree.clade[0].depths().values())
            deep_dist_1 = max(tree.clade[1].depths().values())
            self.assertAlmostEqual(deep_dist_0, deep_dist_1)
        # The midpoint can lie on a branch adjacent to the root; splitting
        # that branch must not inflate the tree or any pairwise distance.
        tree = Phylo.read(StringIO("(A:1.0,B:2.0);"), "newick")
        tree.root_at_midpoint()
        self.assertAlmostEqual(tree.total_branch_length(), 3.0)
        self.assertAlmostEqual(tree.distance("A", "B"), 3.0)
        self.assertEqual(len(tree.root.clades), 2)
        lengths = sorted(c.branch_length for c in tree.root.clades)
        self.assertAlmostEqual(lengths[0], 1.5)
        self.assertAlmostEqual(lengths[1], 1.5)
        tree = Phylo.read(StringIO("(t3:0,(t0:0,t2:1):2);"), "newick")
        tree.root_at_midpoint()
        self.assertAlmostEqual(tree.total_branch_length(), 3.0)
        self.assertAlmostEqual(tree.distance("t0", "t3"), 2.0)
        self.assertAlmostEqual(tree.distance("t0", "t2"), 1.0)
        self.assertAlmostEqual(tree.distance("t2", "t3"), 3.0)

    # Magic method
    def test_str(self):
        """Tree.__str__: pretty-print to a string.

        NB: The exact line counts are liable to change if the object
        constructors change.
        """
        for source, count in zip((EX_APAF, EX_BCL2, EX_DIST), (386, 747, 15)):
            tree = Phylo.read(source, "phyloxml")
            output = str(tree)
            self.assertEqual(len(output.splitlines()), count)


class RerootDistanceTests(unittest.TestCase):
    """Rerooting must not change any distance between two tips.

    Rerooting only moves the root, so the total branch length and the path
    length between every pair of tips must be what they were before.
    """

    def random_tree(self, rng, n_tips):
        """Return a random unrooted-style Newick tree with n_tips tips.

        Clades are joined at random, and the last three are put under a
        trifurcating root half of the time, as tree builders often do.
        """
        clades = [
            Newick.Clade(name=f"t{i}", branch_length=rng.choice([0.0, rng.random()]))
            for i in range(n_tips)
        ]
        stop = rng.choice([2, 3])
        while len(clades) > stop:
            rng.shuffle(clades)
            clades.append(
                Newick.Clade(
                    clades=[clades.pop(), clades.pop()],
                    branch_length=rng.choice([0.0, rng.random()]),
                )
            )
        return Newick.Tree(root=Newick.Clade(clades=clades), rooted=False)

    def tip_distances(self, tree):
        """Return the distance between each pair of tips, keyed by name."""
        tips = sorted(tree.get_terminals(), key=lambda tip: tip.name)
        return {
            (a.name, b.name): tree.distance(a, b)
            for i, a in enumerate(tips)
            for b in tips[i + 1 :]
        }

    def assertDistancesKept(self, tree, total, distances, msg):
        """Check the total branch length and every pairwise tip distance."""
        self.assertAlmostEqual(tree.total_branch_length(), total, msg=msg)
        for pair, value in self.tip_distances(tree).items():
            self.assertAlmostEqual(value, distances[pair], msg=f"{msg} {pair}")

    def test_random_root_with_outgroup(self):
        """Tree.root_with_outgroup keeps tip distances on random trees."""
        rng = random.Random(5390)
        for _ in range(100):
            tree = self.random_tree(rng, rng.randint(4, 10))
            total = tree.total_branch_length()
            distances = self.tip_distances(tree)
            newick = tree.format("newick").strip()
            # Reroot repeatedly on a random clade, putting the new root at a
            # random point along the branch above it or giving no length.
            for _ in range(10):
                clade = rng.choice(list(tree.find_clades()))
                blen = clade.branch_length or 0.0
                tree.root_with_outgroup(
                    clade,
                    outgroup_branch_length=rng.choice([None, rng.random() * blen]),
                )
                self.assertDistancesKept(tree, total, distances, newick)

    def test_random_root_at_midpoint(self):
        """Tree.root_at_midpoint keeps tip distances on random trees."""
        rng = random.Random(5390)
        for _ in range(300):
            tree = self.random_tree(rng, rng.randint(4, 10))
            total = tree.total_branch_length()
            distances = self.tip_distances(tree)
            longest = max(distances.values())
            if not longest:
                # All tips coincide, so there is no midpoint to find
                continue
            newick = tree.format("newick").strip()
            tree.root_at_midpoint()
            self.assertDistancesKept(tree, total, distances, newick)
            # The root is at the midpoint of the longest path between tips
            self.assertEqual(len(tree.root.clades), 2, msg=newick)
            for child in tree.root.clades:
                # depths() counts the child's own branch from the root
                depth = max(child.depths().values())
                self.assertAlmostEqual(depth, longest / 2, msg=newick)

    def test_distance_tree_root_at_midpoint(self):
        """Midpoint rooting an NJ or UPGMA tree keeps its tip distances."""
        rng = random.Random(5390)
        names = [f"s{i}" for i in range(6)]
        constructor = DistanceTreeConstructor()
        for method in (constructor.nj, constructor.upgma):
            for _ in range(60):
                points = [[rng.random() for _ in range(3)] for _ in names]
                matrix = [
                    [math.dist(points[i], points[j]) for j in range(i + 1)]
                    for i in range(len(names))
                ]
                tree = method(DistanceMatrix(names, matrix))
                total = tree.total_branch_length()
                distances = self.tip_distances(tree)
                newick = tree.format("newick").strip()
                tree.root_at_midpoint()
                self.assertDistancesKept(tree, total, distances, newick)


def newick_support(clade):
    """Return the support of a Newick clade."""
    return clade.confidence


def phyloxml_support(clade):
    """Return the confidences of a PhyloXML clade as (value, type) pairs."""
    return tuple((float(conf), conf.type) for conf in clade.confidences) or None


class RerootSupportTests(unittest.TestCase):
    """Rerooting moves support values with their branches.

    A support value describes the split (bipartition) made by the branch
    above a clade, so root_with_outgroup must move it with that branch, as it
    moves the branch length. See upstream issue #745.
    """

    def split_supports(self, tree, get_support=newick_support):
        """Map each split of the tree to the support of its branches.

        A split is named by the taxa on the side of its branch away from the
        first taxon by name. Both branches of a bifurcating root make the same
        split, so where both carry a support value they must agree.
        """
        taxa = frozenset(leaf.name for leaf in tree.get_terminals())
        first = min(taxa)
        supports = {}
        for clade in tree.find_clades():
            if clade is tree.root:
                continue
            side = frozenset(leaf.name for leaf in clade.get_terminals())
            if first in side:
                side = taxa - side
            value = get_support(clade)
            if value is not None:
                self.assertEqual(supports.setdefault(side, value), value, side)
        return supports

    def assert_reroots_keep_splits(self, tree, get_support=newick_support):
        """Reroot at every clade in turn, checking each split keeps its support."""
        before = self.split_supports(tree, get_support)
        self.assertTrue(before)
        for outgroup_branch_length in (None, 0.25):
            for clade in list(tree.find_clades()):
                tree.root_with_outgroup(
                    clade, outgroup_branch_length=outgroup_branch_length
                )
                self.assertEqual(self.split_supports(tree, get_support), before)
                self.assertIsNone(get_support(tree.root))

    def test_report_tree(self):
        """Support moves off the outgroup's path onto the reversed branches."""
        # The second tree from upstream issue #745
        tree = Phylo.read(StringIO("(((A,B)90:1,(C,D)80:1):1,E:1);"), "newick")
        before = self.split_supports(tree)
        self.assertEqual(before, {frozenset("CDE"): 90, frozenset("CD"): 80})
        tree.root_with_outgroup("A")
        self.assertEqual(self.split_supports(tree), before)
        # The 90 used to stay on (A,B), whose branch now splits A from the rest
        self.assertIsNone(tree.common_ancestor("B", "C").confidence)
        self.assertEqual(tree.common_ancestor("C", "E").confidence, 90)
        self.assertIsNone(tree.root.confidence)

    def test_published_test_tree(self):
        """Rerooting the test tree of Czech et al. 2017 on X labels every split."""
        # From https://github.com/jhcepas/test_branch_support_after_tree_rerooting
        tree = Phylo.read(
            StringIO(
                "((C:1.0,D:1.0)10:0.001,(A:1.0,(B:1.0,X:1.0)30:0.01)20:0.1,E:1.0)0:0.0;"
            ),
            "newick",
        )
        before = self.split_supports(tree)
        self.assertEqual(
            before, {frozenset("CD"): 10, frozenset("CDE"): 20, frozenset("BX"): 30}
        )
        tree.root_with_outgroup("X")
        self.assertEqual(self.split_supports(tree), before)
        # The root's own label, 0, used to land on the {C,D,E} branch
        self.assertEqual(tree.common_ancestor("C", "E").confidence, 20)
        self.assertIsNone(tree.root.confidence)

    def test_old_root_branches_merge(self):
        """The old bifurcating root's two branches merge into one."""
        # The first tree from upstream issue #745 gives the split {A,B}|{C,D}
        # two values. The clade below the merged branch keeps its own value,
        # whichever side the outgroup is on, and even if its value is 0.
        for newick, outgroup, ingroup, kept in (
            ("((A,B)90:1,(C,D)80:1);", "A", ("C", "D"), 80),
            ("((A,B)90:1,(C,D)80:1);", "C", ("A", "B"), 90),
            ("((A,B)90:1,(C,D)0:1);", "A", ("C", "D"), 0),
        ):
            tree = Phylo.read(StringIO(newick), "newick")
            with self.assertWarns(BiopythonWarning):
                tree.root_with_outgroup(outgroup)
            self.assertEqual(tree.common_ancestor(*ingroup).confidence, kept)
            self.assertEqual(len(self.split_supports(tree)), 1)
            # The root's branches split the outgroup from the rest: no support
            for clade in tree.root:
                self.assertIsNone(clade.confidence)
        # With a value on only one of the two branches, it is kept, even 0
        for value in (90, 0):
            tree = Phylo.read(StringIO(f"((A,B){value}:1,(C,D):1);"), "newick")
            tree.root_with_outgroup("A")
            self.assertEqual(tree.common_ancestor("C", "D").confidence, value)
            self.assertIsNone(tree.common_ancestor("B", "C").confidence)

    def assertWarnsDropped(self, kept, dropped, outgroup, **kwargs):
        """Check rerooting warns that one of the root's two supports is dropped.

        The outgroup is a tuple of taxon names.
        """
        with self.assertWarns(BiopythonWarning) as cm:
            self.tree.root_with_outgroup(*outgroup, **kwargs)
        message = str(cm.warning)
        self.assertIn(f"keeping {kept} and dropping {dropped}.", message)
        # Reported at the caller's line, not inside Bio.Phylo
        self.assertEqual(cm.filename, __file__)

    def assertNoWarning(self, outgroup, **kwargs):
        """Check rerooting with a tuple of taxon names emits no warning."""
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            self.tree.root_with_outgroup(*outgroup, **kwargs)

    def test_old_root_conflict_warns(self):
        """Dropping one of two different supports on the old root warns."""
        # The ingroup keeps its own value and the outgroup side's is dropped,
        # whichever side the outgroup is on, with or without a branch length,
        # and for an internal outgroup that becomes the root. A 0 is a value.
        for newick, outgroup, kwargs, kept, dropped in (
            ("((A,B)90:1,(C,D)80:1);", ("A",), {}, 80, 90),
            ("((A,B)90:1,(C,D)80:1);", ("C",), {}, 90, 80),
            ("((A,B)90:1,(C,D)80:1);", ("A", "B"), {}, 80, 90),
            (
                "((A,B)90:1,(C,D)80:1);",
                ("A", "B"),
                {"outgroup_branch_length": 0.5},
                80,
                90,
            ),
            ("((A,B)90:1,(C,D)0:1);", ("A",), {}, 0, 90),
            ("((A,B)0:1,(C,D)80:1);", ("A",), {}, 80, 0),
        ):
            self.tree = Phylo.read(StringIO(newick), "newick")
            self.assertWarnsDropped(kept, dropped, outgroup, **kwargs)
        # No warning when the values agree or only one branch has one
        for newick in (
            "((A,B)90:1,(C,D)90:1);",
            "((A,B)90:1,(C,D):1);",
            "((A,B):1,(C,D)90:1);",
            "((A,B)0:1,(C,D):1);",
            "((A,B):1,(C,D):1);",
        ):
            for outgroup, kwargs in (
                (("A",), {}),
                (("C",), {}),
                (("A", "B"), {}),
                (("A", "B"), {"outgroup_branch_length": 0.5}),
            ):
                self.tree = Phylo.read(StringIO(newick), "newick")
                self.assertNoWarning(outgroup, **kwargs)
        # A polytomous old root stays in the tree, so nothing merges
        self.tree = Phylo.read(StringIO("((A,B)90:1,(C,D)80:1,E:1);"), "newick")
        self.assertNoWarning(("A",))
        # Rerooting again only merges the halves of one branch, which agree
        self.tree = Phylo.read(StringIO("((A,B)90:1,(C,D)80:1);"), "newick")
        self.assertWarnsDropped(80, 90, ("A",))
        for outgroup in (("C",), ("B",), ("A", "B"), ("C", "D")):
            self.assertNoWarning(outgroup)
            self.assertNoWarning(outgroup, outgroup_branch_length=0.5)
        # NaN on both branches drops nothing that differs
        for left, right, warns in ((math.nan, math.nan, False), (math.nan, 0.9, True)):
            self.tree = Phylo.BaseTree.Tree(
                Phylo.BaseTree.Clade(
                    clades=[
                        Phylo.BaseTree.Clade(
                            branch_length=1,
                            confidence=left,
                            clades=[
                                Phylo.BaseTree.Clade(name="A"),
                                Phylo.BaseTree.Clade(name="B"),
                            ],
                        ),
                        Phylo.BaseTree.Clade(
                            branch_length=1,
                            confidence=right,
                            clades=[
                                Phylo.BaseTree.Clade(name="C"),
                                Phylo.BaseTree.Clade(name="D"),
                            ],
                        ),
                    ]
                )
            )
            if warns:
                self.assertWarnsDropped(right, left, ("A",))
            else:
                self.assertNoWarning(("A",))

    def test_old_root_conflict_warns_midpoint(self):
        """root_at_midpoint warns once, at the caller's line."""
        tree = Phylo.read(StringIO("((A:1,B:5)90:1,(C:1,D:1)80:1);"), "newick")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            tree.root_at_midpoint()
        self.assertEqual(len(caught), 1)
        self.assertIs(caught[0].category, BiopythonWarning)
        self.assertIn("keeping", str(caught[0].message))
        self.assertEqual(caught[0].filename, __file__)
        # Equal values do not warn
        tree = Phylo.read(StringIO("((A:1,B:5)90:1,(C:1,D:1)90:1);"), "newick")
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            tree.root_at_midpoint()

    def test_old_root_conflict_warns_phyloxml(self):
        """PhyloXML supports differ if their values or their types differ."""
        xml = """\
<phyloxml xmlns="http://www.phyloxml.org">
  <phylogeny rooted="true">
    <clade>
      <clade>
        <branch_length>1</branch_length>
        <confidence type="bootstrap">90</confidence>
        <clade><name>A</name></clade>
        <clade><name>B</name></clade>
      </clade>
      <clade>
        <branch_length>1</branch_length>
        %s
        <clade><name>C</name></clade>
        <clade><name>D</name></clade>
      </clade>
    </clade>
  </phylogeny>
</phyloxml>
"""
        for confidences in ('<confidence type="bootstrap">90</confidence>', ""):
            self.tree = Phylo.read(StringIO(xml % confidences), "phyloxml")
            self.assertNoWarning(("A",))
        # The ingroup (C,D) keeps its own list; the message gives the types
        for confidences, kept in (
            ('<confidence type="bootstrap">80</confidence>', "[bootstrap 80.0]"),
            ('<confidence type="probability">90</confidence>', "[probability 90.0]"),
            (
                '<confidence type="bootstrap">90</confidence>'
                '<confidence type="probability">0.95</confidence>',
                "[bootstrap 90.0, probability 0.95]",
            ),
        ):
            self.tree = Phylo.read(StringIO(xml % confidences), "phyloxml")
            self.assertWarnsDropped(kept, "[bootstrap 90.0]", ("A",))
        # The order of the confidence elements carries no meaning
        two = xml.replace(
            '<confidence type="bootstrap">90</confidence>',
            '<confidence type="bootstrap">90</confidence>'
            '<confidence type="probability">0.95</confidence>',
        )
        self.tree = Phylo.read(
            StringIO(
                two
                % (
                    '<confidence type="probability">0.95</confidence>'
                    '<confidence type="bootstrap">90</confidence>'
                )
            ),
            "phyloxml",
        )
        self.assertNoWarning(("A",))

    def test_internal_outgroup_becomes_root(self):
        """An internal outgroup's support moves to the branch below it."""
        tree = Phylo.read(StringIO("((A,B)90:1,(C,D)80:1,E:1)100;"), "newick")
        before = self.split_supports(tree)
        tree.root_with_outgroup("A", "B")
        self.assertIs(tree.root, tree.common_ancestor("A", "B"))
        self.assertIsNone(tree.root.confidence)
        self.assertEqual(tree.common_ancestor("C", "E").confidence, 90)
        self.assertEqual(self.split_supports(tree), before)

    def test_outgroup_branch_length(self):
        """Both halves of the branch holding the new root carry its support."""
        tree = Phylo.read(StringIO("(((A,B)90:1,(C,D)80:1)70:1,E:1,F:1);"), "newick")
        before = self.split_supports(tree)
        tree.root_with_outgroup("A", "B", outgroup_branch_length=0.5)
        self.assertEqual([clade.confidence for clade in tree.root], [90, 90])
        self.assertEqual(tree.common_ancestor("E", "F").confidence, 70)
        self.assertIsNone(tree.root.confidence)
        self.assertEqual(self.split_supports(tree), before)

    def test_direct_child_outgroup_branch_length(self):
        """A new root on the old root's merged branch puts one value on both halves."""
        # Both of the old root's branches make the split {A,B}|{C,D}. The
        # ingroup keeps its own value, and the outgroup's half gets a copy.
        for outgroup, kept in ((("A", "B"), 80), (("C", "D"), 90)):
            tree = Phylo.read(StringIO("((A:1,B:1)90:1,(C:1,D:1)80:1);"), "newick")
            with self.assertWarns(BiopythonWarning):
                tree.root_with_outgroup(*outgroup, outgroup_branch_length=0.5)
            self.assertEqual([clade.confidence for clade in tree.root], [kept, kept])
            self.assertIsNone(tree.root.confidence)
        # With a value on only one of the two branches, both halves get it
        for newick in ("((A:1,B:1)90:1,(C:1,D:1):1);", "((A:1,B:1):1,(C:1,D:1)90:1);"):
            tree = Phylo.read(StringIO(newick), "newick")
            tree.root_with_outgroup("A", "B", outgroup_branch_length=0.5)
            self.assertEqual([clade.confidence for clade in tree.root], [90, 90])
        # PhyloXML: the copies are separate lists of separate objects
        xml = """\
<phyloxml xmlns="http://www.phyloxml.org">
  <phylogeny rooted="true">
    <clade>
      <clade>
        <branch_length>1</branch_length>
        <confidence type="bootstrap">90</confidence>
        <clade><name>A</name></clade>
        <clade><name>B</name></clade>
      </clade>
      <clade>
        <branch_length>1</branch_length>
        <confidence type="bootstrap">80</confidence>
        <clade><name>C</name></clade>
        <clade><name>D</name></clade>
      </clade>
    </clade>
  </phylogeny>
</phyloxml>
"""
        tree = Phylo.read(StringIO(xml), "phyloxml")
        with self.assertWarns(BiopythonWarning):
            tree.root_with_outgroup("A", "B", outgroup_branch_length=0.5)
        first, second = tree.root.clades
        self.assertEqual(phyloxml_support(first), ((80, "bootstrap"),))
        self.assertEqual(phyloxml_support(second), ((80, "bootstrap"),))
        self.assertIsNot(first.confidences, second.confidences)
        self.assertIsNot(first.confidences[0], second.confidences[0])

    def test_root_at_midpoint(self):
        """Tree.root_at_midpoint moves support values too."""
        tree = Phylo.read(StringIO("((A:1,B:1)90:4,(C:1,D:1)80:1,E:1);"), "newick")
        before = self.split_supports(tree)
        tree.root_at_midpoint()
        # The midpoint is on the branch of (A,B)
        self.assertEqual([clade.confidence for clade in tree.root], [90, 90])
        self.assertEqual(tree.common_ancestor("C", "D").confidence, 80)
        self.assertIsNone(tree.root.confidence)
        self.assertEqual(self.split_supports(tree), before)

    def test_reroot_everywhere(self):
        """Each split keeps its support through repeated rerooting."""
        for newick in (
            "((A,B)1,((C,D)2,(E,F)3)4,(G,H)5);",
            "(((A,B)1,((C,D)2,(E,F)3)4)6,(G,H));",
            "((A,B,C),(D,(E,F)2)3);",
        ):
            for rooted in (False, True):
                tree = Phylo.read(StringIO(newick), "newick", rooted=rooted)
                self.assert_reroots_keep_splits(tree)
                self.assertTrue(tree.rooted)

    def test_phyloxml_confidences(self):
        """A PhyloXML clade's whole list of confidences moves with its branch."""
        xml = """\
<phyloxml xmlns="http://www.phyloxml.org">
  <phylogeny rooted="false">
    <clade>
      <clade>
        <branch_length>1</branch_length>
        <confidence type="bootstrap">90</confidence>
        <confidence type="probability">0.95</confidence>
        <clade><name>A</name></clade>
        <clade><name>B</name></clade>
      </clade>
      <clade>
        <branch_length>1</branch_length>
        <confidence type="bootstrap">80</confidence>
        <clade><name>C</name></clade>
        <clade><name>D</name></clade>
      </clade>
      <clade><name>E</name></clade>
    </clade>
  </phylogeny>
</phyloxml>
"""
        tree = Phylo.read(StringIO(xml), "phyloxml")
        tree.root_with_outgroup("A")
        moved = tree.common_ancestor("C", "E")
        self.assertEqual(
            phyloxml_support(moved), ((90, "bootstrap"), (0.95, "probability"))
        )
        self.assertEqual(tree.common_ancestor("B", "C").confidences, [])
        self.assertTrue(tree.rooted)
        # Both halves of the outgroup's branch get the values, in separate lists
        tree = Phylo.read(StringIO(xml), "phyloxml")
        tree.root_with_outgroup("A", "B", outgroup_branch_length=0.5)
        first, second = tree.root.clades
        self.assertEqual(phyloxml_support(first), phyloxml_support(second))
        self.assertEqual(len(first.confidences), 2)
        self.assertIsNot(first.confidences, second.confidences)
        self.assertIsNot(first.confidences[0], second.confidences[0])
        self.assertEqual(tree.root.confidences, [])

    def test_phyloxml_file(self):
        """Each split of a PhyloXML file keeps its support through rerooting."""
        tree = Phylo.read(EX_APAF, "phyloxml")
        self.assert_reroots_keep_splits(tree, phyloxml_support)


class MixinTests(unittest.TestCase):
    """Tests for TreeMixin methods."""

    def setUp(self):
        self.phylogenies = list(Phylo.parse(EX_PHYLO, "phyloxml"))

    # Traversal methods

    def test_find_elements(self):
        """TreeMixin: find_elements() method."""
        # From the docstring example
        tree = self.phylogenies[5]
        matches = list(tree.find_elements(PhyloXML.Taxonomy, code="OCTVU"))
        self.assertEqual(len(matches), 1)
        self.assertIsInstance(matches[0], PhyloXML.Taxonomy)
        self.assertEqual(matches[0].code, "OCTVU")
        self.assertEqual(matches[0].scientific_name, "Octopus vulgaris")
        # Iteration and regexps
        tree = self.phylogenies[10]
        for point, alt in zip(
            tree.find_elements(geodetic_datum=r"WGS\d{2}"), (472, 10, 452)
        ):
            self.assertIsInstance(point, PhyloXML.Point)
            self.assertEqual(point.geodetic_datum, "WGS84")
            self.assertAlmostEqual(point.alt, alt)
        # class filter
        tree = self.phylogenies[4]
        events = list(tree.find_elements(PhyloXML.Events))
        self.assertEqual(len(events), 2)
        self.assertEqual(events[0].speciations, 1)
        self.assertEqual(events[1].duplications, 1)
        # string filter & find_any
        tree = self.phylogenies[3]
        taxonomy = tree.find_any("B. subtilis")
        self.assertEqual(taxonomy.scientific_name, "B. subtilis")
        # integer filter
        tree = Phylo.read(EX_APAF, "phyloxml")
        domains = list(tree.find_elements(start=5))
        self.assertEqual(len(domains), 8)
        for dom in domains:
            self.assertEqual(dom.start, 5)
            self.assertEqual(dom.value, "CARD")

    def test_find_clades(self):
        """TreeMixin: find_clades() method."""
        # boolean filter
        for clade, name in zip(
            self.phylogenies[10].find_clades(name=True), list("ABCD")
        ):
            self.assertIsInstance(clade, PhyloXML.Clade)
            self.assertEqual(clade.name, name)
        # finding deeper attributes
        octo = list(self.phylogenies[5].find_clades(code="OCTVU"))
        self.assertEqual(len(octo), 1)
        self.assertIsInstance(octo[0], PhyloXML.Clade)
        self.assertEqual(octo[0].taxonomies[0].code, "OCTVU")
        # string filter
        dee = next(self.phylogenies[10].find_clades("D"))
        self.assertEqual(dee.name, "D")

    def test_find_terminal(self):
        """TreeMixin: find_elements() with terminal argument."""
        for tree, total, extern, intern in zip(
            self.phylogenies,
            (6, 6, 7, 18, 21, 27, 7, 9, 9, 19, 15, 9, 6),
            (3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 3, 3),
            (3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
        ):
            self.assertEqual(len(list(tree.find_elements())), total)
            self.assertEqual(len(list(tree.find_elements(terminal=True))), extern)
            self.assertEqual(len(list(tree.find_elements(terminal=False))), intern)

    def test_get_path(self):
        """TreeMixin: get_path() method."""
        path = self.phylogenies[1].get_path("B")
        self.assertEqual(len(path), 2)
        self.assertAlmostEqual(path[0].branch_length, 0.06)
        self.assertAlmostEqual(path[1].branch_length, 0.23)
        self.assertEqual(path[1].name, "B")

    def test_trace(self):
        """TreeMixin: trace() method."""
        tree = self.phylogenies[1]
        path = tree.trace("A", "C")
        self.assertEqual(len(path), 3)
        self.assertAlmostEqual(path[0].branch_length, 0.06)
        self.assertAlmostEqual(path[2].branch_length, 0.4)
        self.assertEqual(path[2].name, "C")

    def test_missing_target_trace(self):
        """TreeMixin: trace() method with missing target."""
        tree = self.phylogenies[1]
        with self.assertRaises(ValueError) as cm:
            tree.trace("Aaa", "C")
        self.assertEqual(str(cm.exception), "target 'Aaa' is not in this tree")

    # Information methods

    def test_common_ancestor(self):
        """TreeMixin: common_ancestor() method."""
        tree = self.phylogenies[1]
        lca = tree.common_ancestor("A", "B")
        self.assertEqual(lca, tree.clade[0])
        lca = tree.common_ancestor("A", "C")
        self.assertEqual(lca, tree.clade)
        tree = self.phylogenies[10]
        lca = tree.common_ancestor("A", "B", "C")
        self.assertEqual(lca, tree.clade[0])

    def test_depths(self):
        """TreeMixin: depths() method."""
        tree = self.phylogenies[1]
        depths = tree.depths()
        self.assertEqual(len(depths), 5)
        for found, expect in zip(
            sorted(depths.values()), [0, 0.060, 0.162, 0.290, 0.400]
        ):
            self.assertAlmostEqual(found, expect)

    def test_distance(self):
        """TreeMixin: distance() method."""
        t = self.phylogenies[1]
        self.assertAlmostEqual(t.distance("A"), 0.162)
        self.assertAlmostEqual(t.distance("B"), 0.29)
        self.assertAlmostEqual(t.distance("C"), 0.4)
        self.assertAlmostEqual(t.distance("A", "B"), 0.332)
        self.assertAlmostEqual(t.distance("A", "C"), 0.562)
        self.assertAlmostEqual(t.distance("B", "C"), 0.69)

    def test_is_bifurcating(self):
        """TreeMixin: is_bifurcating() method."""
        for tree, is_b in zip(
            self.phylogenies, (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1)
        ):
            self.assertEqual(tree.is_bifurcating(), is_b)

    def test_is_monophyletic(self):
        """TreeMixin: is_monophyletic() method."""
        tree = self.phylogenies[10]
        abcd = tree.get_terminals()
        abc = tree.clade[0].get_terminals()
        ab = abc[:2]
        d = tree.clade[1].get_terminals()
        self.assertEqual(tree.is_monophyletic(abcd), tree.root)
        self.assertEqual(tree.is_monophyletic(abc), tree.clade[0])
        self.assertFalse(tree.is_monophyletic(ab))
        self.assertEqual(tree.is_monophyletic(d), tree.clade[1])
        # Alternate argument form
        self.assertEqual(tree.is_monophyletic(*abcd), tree.root)

    def test_total_branch_length(self):
        """TreeMixin: total_branch_length() method."""
        tree = self.phylogenies[1]
        self.assertAlmostEqual(tree.total_branch_length(), 0.792)
        self.assertAlmostEqual(tree.clade[0].total_branch_length(), 0.392)

    # Tree manipulation methods

    def test_collapse(self):
        """TreeMixin: collapse() method."""
        tree = self.phylogenies[1]
        parent = tree.collapse(tree.clade[0])
        self.assertEqual(len(parent), 3)
        for clade, name, blength in zip(parent, ("C", "A", "B"), (0.4, 0.162, 0.29)):
            self.assertEqual(clade.name, name)
            self.assertAlmostEqual(clade.branch_length, blength)

    def test_collapse_all(self):
        """TreeMixin: collapse_all() method."""
        tree = Phylo.read(EX_APAF, "phyloxml")
        d1 = tree.depths()
        tree.collapse_all()
        d2 = tree.depths()
        # Total branch lengths should not change
        for clade in d2:
            self.assertAlmostEqual(d1[clade], d2[clade])
        # No internal nodes should remain except the root
        self.assertEqual(len(tree.get_terminals()), len(tree.clade))
        self.assertEqual(len(list(tree.find_clades(terminal=False))), 1)
        # Again, with a target specification
        tree = Phylo.read(EX_APAF, "phyloxml")
        d1 = tree.depths()
        internal_node_ct = len(tree.get_nonterminals())
        tree.collapse_all(lambda c: c.branch_length < 0.1)  # noqa: E731
        d2 = tree.depths()
        # Should have collapsed 7 internal nodes
        self.assertEqual(len(tree.get_nonterminals()), internal_node_ct - 7)
        for clade in d2:
            self.assertAlmostEqual(d1[clade], d2[clade])

    def test_collapse_without_branch_lengths(self):
        tree = Phylo.read(StringIO("(((A,B),C),D);"), "newick")
        tree.collapse_all()
        self.assertEqual(len(tree.get_terminals()), 4)
        self.assertEqual(len(tree.get_nonterminals()), 1)

    def test_ladderize(self):
        """TreeMixin: ladderize() method."""

        def ordered_names(tree):
            return [n.name for n in tree.get_terminals()]

        tree = self.phylogenies[10]
        self.assertEqual(ordered_names(tree), list("ABCD"))
        tree.ladderize()
        self.assertEqual(ordered_names(tree), list("DABC"))
        tree.ladderize(reverse=True)
        self.assertEqual(ordered_names(tree), list("ABCD"))

    def test_prune(self):
        """TreeMixin: prune() method."""
        tree = self.phylogenies[10]
        # Taxon in a trifurcation -- no collapse afterward
        parent = tree.prune(name="B")
        self.assertEqual(len(parent.clades), 2)
        self.assertEqual(parent.clades[0].name, "A")
        self.assertEqual(parent.clades[1].name, "C")
        self.assertEqual(len(tree.get_terminals()), 3)
        self.assertEqual(len(tree.get_nonterminals()), 2)
        # Taxon in a bifurcation -- collapse
        tree = self.phylogenies[0]
        parent = tree.prune(name="A")
        self.assertEqual(len(parent.clades), 2)
        for clade, name, blen in zip(parent, "BC", (0.29, 0.4)):
            self.assertTrue(clade.is_terminal())
            self.assertEqual(clade.name, name)
            self.assertAlmostEqual(clade.branch_length, blen)
        self.assertEqual(len(tree.get_terminals()), 2)
        self.assertEqual(len(tree.get_nonterminals()), 1)
        # Taxon just below the root -- don't screw up
        tree = self.phylogenies[1]
        parent = tree.prune(name="C")
        self.assertEqual(parent, tree.root)
        self.assertEqual(tree.root.branch_length, 0.06)
        self.assertEqual(len(parent.clades), 2)
        for clade, name, blen in zip(parent, "AB", (0.102, 0.23)):
            self.assertTrue(clade.is_terminal())
            self.assertEqual(clade.name, name)
            self.assertAlmostEqual(clade.branch_length, blen)
        self.assertEqual(len(tree.get_terminals()), 2)
        self.assertEqual(len(tree.get_nonterminals()), 1)
        # Taxon just below root
        tree = self.phylogenies[13]
        parent = tree.prune(name="C")
        self.assertEqual(tree.root.branch_length, 0.1)

    def test_split(self):
        """TreeMixin: split() method."""
        tree = self.phylogenies[0]
        C = tree.clade[1]
        C.split()
        self.assertEqual(len(C), 2)
        self.assertEqual(len(tree.get_terminals()), 4)
        self.assertEqual(len(tree.get_nonterminals()), 3)
        C[0].split(3, 0.5)
        self.assertEqual(len(tree.get_terminals()), 6)
        self.assertEqual(len(tree.get_nonterminals()), 4)
        for clade, name, blen in zip(C[0], ("C00", "C01", "C02"), (0.5, 0.5, 0.5)):
            self.assertTrue(clade.is_terminal())
            self.assertEqual(clade.name, name)
            self.assertEqual(clade.branch_length, blen)

    def test_resolve_polytomies(self):
        """TreeMixin: resolve_polytomies() method."""
        poly_tree_str = "(A:2,(B:1,C:2,D:4,(E:3,(F:2,G:3,H:1):2):1):5);\n"
        #                         ^   ^             ^   polytomies
        # Each polytomy resolves to a balanced hierarchy of zero-length
        # bifurcations, preserving the child order left to right.
        resolved_tree_str = (
            "(A:2,((B:1,C:2):0,(D:4,(E:3,(F:2,(G:3,H:1):0):2):1):0):5);\n"
        )
        poly_tree = Phylo.read(StringIO(poly_tree_str), "newick")
        resolved_tree = Phylo.read(StringIO(resolved_tree_str), "newick")
        self.assertFalse(poly_tree.is_bifurcating())
        poly_tree.resolve_polytomies(shuffle=False, recursive=True)
        self.assertTrue(poly_tree.is_bifurcating())
        self.assertEqual(poly_tree.format("newick"), resolved_tree.format("newick"))

    def test_resolve_polytomies_wide(self):
        """TreeMixin: resolve_polytomies() on a star wider than the recursion limit."""
        # A star tree, as produced by e.g. FastTree on unresolvable data.
        # 1500 children exceeds Python's default recursion limit of 1000;
        # a chain-per-child resolution either dies while building (the
        # original recursive implementation) or leaves a tree too deep for
        # the library's recursive traversals to validate or serialise.
        # The balanced resolution keeps the tree ~log2(n) deep, so the
        # ordinary machinery must work on the result.
        size = 1500
        star = "(" + ",".join("t%d:1" % i for i in range(size)) + ");"
        tree = Phylo.read(StringIO(star), "newick")
        tree.resolve_polytomies()
        self.assertTrue(tree.is_bifurcating())
        for clade in tree.find_clades():
            self.assertLessEqual(len(clade.clades), 2)
        # Serialisation round trip, preserving leaf count and order
        reread = Phylo.read(StringIO(tree.format("newick")), "newick")
        self.assertEqual(
            [leaf.name for leaf in reread.get_terminals()],
            ["t%d" % i for i in range(size)],
        )


# ---------------------------------------------------------

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
