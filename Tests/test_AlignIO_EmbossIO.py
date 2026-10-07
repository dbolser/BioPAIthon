# Copyright 2008-2014 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.AlignIO.EmbossIO module."""

import unittest
from io import StringIO

import support

from Bio.AlignIO.EmbossIO import EmbossIterator

# http://emboss.sourceforge.net/docs/themes/alnformats/align.simple
with open(support.DATA / "Emboss" / "alignret.txt") as handle:
    simple_example = handle.read()

# http://emboss.sourceforge.net/docs/themes/alnformats/align.pair
with open(support.DATA / "Emboss" / "water.txt") as handle:
    pair_example = handle.read()

with open(support.DATA / "Emboss" / "needle.txt") as handle:
    pair_example2 = handle.read()

with open(support.DATA / "Emboss" / "needle_overhang.txt") as handle:
    pair_example3 = handle.read()


class TestEmbossIO(unittest.TestCase):
    def test_pair_example(self):
        alignments = list(EmbossIterator(StringIO(pair_example)))
        self.assertEqual(len(alignments), 1)
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual([r.id for r in alignments[0]], ["IXI_234", "IXI_235"])

    def test_simple_example(self):
        alignments = list(EmbossIterator(StringIO(simple_example)))
        self.assertEqual(len(alignments), 1)
        self.assertEqual(len(alignments[0]), 4)
        self.assertEqual(
            [r.id for r in alignments[0]], ["IXI_234", "IXI_235", "IXI_236", "IXI_237"]
        )

    def test_pair_plus_simple(self):
        alignments = list(EmbossIterator(StringIO(pair_example + simple_example)))
        self.assertEqual(len(alignments), 2)
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(len(alignments[1]), 4)
        self.assertEqual([r.id for r in alignments[0]], ["IXI_234", "IXI_235"])
        self.assertEqual(
            [r.id for r in alignments[1]], ["IXI_234", "IXI_235", "IXI_236", "IXI_237"]
        )

    def test_pair_example2(self):
        alignments = list(EmbossIterator(StringIO(pair_example2)))
        self.assertEqual(len(alignments), 5)
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(
            [r.id for r in alignments[0]], ["ref_rec", "gi|94968718|receiver"]
        )
        self.assertEqual(
            [r.id for r in alignments[4]], ["ref_rec", "gi|94970041|receiver"]
        )

    def test_pair_example3(self):
        alignments = list(EmbossIterator(StringIO(pair_example3)))
        self.assertEqual(len(alignments), 1)
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual([r.id for r in alignments[0]], ["asis", "asis"])


class TestEmbossIOMalformed(unittest.TestCase):
    """Malformed EMBOSS input must raise ValueError, even under python -O."""

    header = """\
#=======================================
#
# Aligned_sequences: 2
# 1: seqA
# 2: seqB
#
# Length: 20
#
#=======================================

"""
    first_block = """\
seqA               1 ACGTACGTAC     10
                     |||||||||
seqB               1 ACGTACGTA-      9

"""

    def parse(self, text):
        return list(EmbossIterator(StringIO(text)))

    def test_all_gap_line_at_digit_boundary(self):
        """An all gap line may have a start with more digits than its end.

        EMBOSS writes "10 ---------- 9" for a gap-only line after nine
        letters. Comparing the start and end as strings put "10" before "9",
        so this valid line used to fail an assert.
        """
        text = (
            self.header
            + self.first_block
            + "seqA              11 ACGTACGTAC     20\n"
            + "\n"
            + "seqB              10 ----------      9\n"
        )
        alignments = self.parse(text)
        self.assertEqual(len(alignments), 1)
        self.assertEqual(alignments[0][1].seq, "ACGTACGTA-----------")

    def test_wrong_identifier(self):
        """Each sequence line must carry the expected (maybe truncated) id."""
        text = (
            self.header
            + self.first_block
            + "seqA              11 ACGTACGTAC     20\n"
            + "\n"
            + "seqC              10 A---------     10\n"
        )
        with self.assertRaises(ValueError) as cm:
            self.parse(text)
        self.assertIn("Expected identifier 'seqB'", str(cm.exception))

    def test_identifier_lines_misnumbered(self):
        """The header lists the sequences as 1, 2, ... in order."""
        text = self.header.replace("# 2: seqB", "# 3: seqB") + self.first_block
        with self.assertRaises(ValueError) as cm:
            self.parse(text)
        self.assertIn("Expected identifier line for sequence 2", str(cm.exception))

    def test_identifier_line_not_numbered(self):
        """An identifier line without a number gets the same error."""
        text = self.header.replace("# 2: seqB", "# x: seqB") + self.first_block
        with self.assertRaises(ValueError) as cm:
            self.parse(text)
        self.assertIn("Expected identifier line for sequence 2", str(cm.exception))

    def test_letters_on_line_with_start_after_end(self):
        """Only a line with no letters can have its start after its end."""
        text = (
            self.header
            + self.first_block
            + "seqA              11 ACGTACGTAC     20\n"
            + "\n"
            + "seqB              12 A---------     10\n"
        )
        with self.assertRaises(ValueError) as cm:
            self.parse(text)
        self.assertIn("Expected start before end", str(cm.exception))

    def test_incomplete_last_block(self):
        """The last block needs a line for every sequence."""
        text = (
            self.header + self.first_block + "seqA              11 ACGTACGTAC     20\n"
        )
        with self.assertRaises(ValueError) as cm:
            self.parse(text)
        self.assertIn(
            "Expected a line for each of the 2 sequences in the last block, found 1",
            str(cm.exception),
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
