# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for the Bio.Align global_align and local_align convenience functions.

Upstream Biopython added these tests (#5290) to ``pairwise2_testCases.py``,
each one alongside the ``Bio.pairwise2`` call it mirrors. BioPAIthon has
removed ``Bio.pairwise2``, so only the ``Bio.Align`` halves are kept here;
the comments about ``pairwise2`` defaults explain the parameters chosen.
"""

import unittest

import numpy as np

from Bio import Align
from Bio.Align import substitution_matrices


class TestPairwiseGlobal(unittest.TestCase):
    """Test some usual global alignments."""

    def test_globalxx_simple(self):
        """Test globalxx."""
        # Note that the PairwiseAligner defaults to gap_score = -1.0,
        # while pairwise2 defaults to a zero gap score, so we need to
        # set gap_score explicitly here.
        alignments = Align.global_align("GAACT", "GAT", g=0.0)  # gap_score
        self.assertEqual(len(alignments), 2)
        self.assertAlmostEqual(alignments.score, 3.0)
        self.assertAlmostEqual(alignments[0].score, 3.0)
        self.assertEqual(
            str(alignments[0]),
            """\
target            0 GAACT 5
                  0 ||--| 5
query             0 GA--T 3
""",
        )
        self.assertAlmostEqual(alignments[1].score, 3.0)
        self.assertEqual(
            str(alignments[1]),
            """\
target            0 GAACT 5
                  0 |-|-| 5
query             0 G-A-T 3
""",
        )

    def test_globalxx_simple2(self):
        """Do the same test with sequence order reversed."""
        # Note that the PairwiseAligner defaults to gap_score = -1.0,
        # while pairwise2 defaults to a zero gap score, so we need to
        # set gap_score explicitly here.
        alignments = Align.global_align("GAT", "GAACT", g=0.0)  # gap_score
        self.assertEqual(len(alignments), 2)
        self.assertAlmostEqual(alignments.score, 3.0)
        self.assertAlmostEqual(alignments[0].score, 3.0)
        self.assertEqual(
            str(alignments[0]),
            """\
target            0 GA--T 3
                  0 ||--| 5
query             0 GAACT 5
""",
        )
        self.assertAlmostEqual(alignments[1].score, 3.0)
        self.assertEqual(
            str(alignments[1]),
            """\
target            0 G-A-T 3
                  0 |-|-| 5
query             0 GAACT 5
""",
        )

    def test_list_input(self):
        """Do a global alignment with sequences supplied as lists."""
        # Note that the PairwiseAligner defaults to gap_score = -1.0,
        # while pairwise2 defaults to a zero gap score, so we need to
        # set gap_score explicitly here.
        alignments = Align.global_align(
            ["Gly", "Ala", "Thr"],
            ["Gly", "Ala", "Ala", "Cys", "Thr"],
            g=0.0,  # gap_score
        )
        self.assertEqual(len(alignments), 2)
        self.assertEqual(alignments.score, 3.0)
        alignment = alignments[0]
        self.assertEqual(alignment[0], ["Gly", "Ala", None, None, "Thr"])
        self.assertEqual(alignment[1], ["Gly", "Ala", "Ala", "Cys", "Thr"])
        self.assertEqual(alignment.score, 3.0)
        alignment = alignments[1]
        self.assertEqual(alignment[0], ["Gly", None, "Ala", None, "Thr"])
        self.assertEqual(alignment[1], ["Gly", "Ala", "Ala", "Cys", "Thr"])
        self.assertEqual(alignment.score, 3.0)


class TestPairwiseLocal(unittest.TestCase):
    """Test some simple local alignments."""

    def setUp(self):
        self.blosum62 = Align.substitution_matrices.load("BLOSUM62")

    def test_localxs_1(self):
        """Test localxx."""
        # From Biopython 1.74 on this should only give one alignment, since
        # we disallow leading and trailing 'zero-extensions'
        alignments = Align.local_align("AxBx", "zABz", o=-0.1)  # open_gap_score
        self.assertEqual(len(alignments), 1)
        self.assertEqual(alignments.score, 1.9)
        self.assertEqual(
            str(alignments[0]),
            """\
target            0 AxB 3
                  0 |-| 3
query             1 A-B 3
""",
        )

    def test_localxs_2(self):
        """Test localxx with ``full_sequences=True``."""
        # From Biopython 1.74 on this should only give one alignment, since
        # we disallow leading and trailing 'zero-extensions'
        alignments = Align.local_align("AxBx", "zABz", o=-0.1)  # open_gap_score
        self.assertEqual(len(alignments), 1)
        self.assertEqual(alignments.score, 1.9)
        alignment = alignments[0]
        alignment.coordinates = np.column_stack(([0, 0], alignment.coordinates, [4, 4]))
        self.assertEqual(
            str(alignment),
            """\
target            0 -AxBx 4
                  0 -|-|. 5
query             0 zA-Bz 4
""",
        )

    def test_localds_zero_score_segments_symmetric(self):
        """Test if alignment is independent on direction of sequence."""
        alignments1 = Align.local_align(
            "CWHISLKM",
            "CWHGISGLKM",
            s=self.blosum62,  # substitution_matrix
            o=-11,  # open_gap_score
            x=-1,  # extend_gap_score
        )
        alignments2 = Align.local_align(
            "MKLSIHWC",
            "MKLGSIGHWC",
            s=self.blosum62,  # substitution_matrix
            o=-11,  # open_gap_score
            x=-1,  # extend_gap_score
        )
        # pairwise2.align.localds gave one alignment in each direction.
        self.assertEqual(len(alignments1), 1)
        self.assertEqual(len(alignments2), 1)

    def test_localxs_generic(self):
        """Test the generic method with local alignments."""
        # From Biopython 1.74 on this should only give one alignment, since
        # we disallow leading and trailing 'zero-extensions'
        alignments = Align.local_align("AxBx", "zABz", o=-0.1)  # open_gap_score
        self.assertEqual(len(alignments), 1)
        self.assertEqual(alignments.score, 1.9)
        self.assertEqual(
            str(alignments[0]),
            """\
target            0 AxB 3
                  0 |-| 3
query             1 A-B 3
""",
        )

    def test_localms(self):
        """Two different local alignments."""
        alignments = Align.local_align(
            "xxxABCDxxx",
            "zzzABzzCDz",
            m=(1.0, -0.5),  # match_score, mismatch_score
            o=-3,  # open_gap_score
            x=-1,  # extend_gap_score
        )
        self.assertEqual(len(alignments), 2)
        self.assertEqual(alignments.score, 2.0)
        self.assertEqual(
            str(alignments[0]),
            """\
target            3 AB 5
                  0 || 2
query             3 AB 5
""",
        )
        self.assertEqual(
            str(alignments[1]),
            """\
target            5 CD 7
                  0 || 2
query             7 CD 9
""",
        )

    def test_blosum62(self):
        """Test localds with blosum62."""
        self.assertEqual(1, self.blosum62[("K", "Q")])
        self.assertEqual(4, self.blosum62[("A", "A")])
        self.assertEqual(8, self.blosum62[("H", "H")])
        alignments = Align.local_align(
            "VKAHGKKV",
            "FQAHCAGV",
            s=self.blosum62,  # substitution_matrix
            o=-4,  # open_gap_score
            x=-4,  # extend_gap_score
        )
        self.assertEqual(len(alignments), 1)
        self.assertAlmostEqual(alignments.score, 13)
        self.assertEqual(
            str(alignments[0]),
            """\
target            1 KAH 4
                  0 .|| 3
query             1 QAH 4
""",
        )

    def test_empty_result(self):
        """Return no alignment."""
        # Note that the PairwiseAligner defaults to gap_score = -1.0,
        # while pairwise2 defaults to a zero gap score, so we need to
        # set gap_score explicitly here.
        alignments = Align.local_align("AT", "GC", g=0)  # gap_score
        self.assertEqual(len(alignments), 0)


class TestScoreOnly(unittest.TestCase):
    """Test parameter ``score_only``."""

    def test_score_only_global(self):
        """Test ``score_only`` in a global alignment."""
        # Note that the PairwiseAligner defaults to gap_score = -1.0,
        # while pairwise2 defaults to a zero gap score, so we need to
        # set gap_score explicitly here.
        alignments = Align.global_align("GAACT", "GAT", g=0)  # gap_score
        self.assertAlmostEqual(alignments.score, 3.0)

    def test_score_only_local(self):
        """Test ``score_only`` in a local alignment."""
        alignments = Align.local_align(
            "xxxABCDxxx",
            "zzzABzzCDz",
            m=(1.0, -0.5),  # match_score, mismatch_score
            o=-3,  # open_gap_score
            x=-1,  # extend_gap_score
        )
        self.assertAlmostEqual(alignments.score, 2.0)


class TestPairwiseOpenPenalty(unittest.TestCase):
    """Alignments with gap-open penalty."""

    def test_match_score_open_penalty1(self):
        """Test 1."""
        alignments = Align.global_align(
            "AA",
            "A",
            m=(2.0, -1.0),  # match_score, mismatch_score
            o=-0.1,  # open_gap_score
            x=0,  # extend_gap_score
        )
        self.assertEqual(len(alignments), 2)
        self.assertAlmostEqual(alignments.score, 1.9)
        self.assertEqual(
            str(alignments[0]),
            """\
target            0 AA 2
                  0 -| 2
query             0 -A 1
""",
        )
        self.assertEqual(
            str(alignments[1]),
            """\
target            0 AA 2
                  0 |- 2
query             0 A- 1
""",
        )

    def test_match_score_open_penalty2(self):
        """Test 2."""
        alignments = Align.global_align(
            "GAA",
            "GA",
            m=(1.5, 0),  # match_score, mismatch_score
            o=-0.1,  # open_gap_score
            x=0,  # extend_gap_score
        )
        self.assertEqual(len(alignments), 2)
        self.assertAlmostEqual(alignments.score, 2.9)
        self.assertEqual(
            str(alignments[0]),
            """\
target            0 GAA 3
                  0 |-| 3
query             0 G-A 2
""",
        )
        self.assertEqual(
            str(alignments[1]),
            """\
target            0 GAA 3
                  0 ||- 3
query             0 GA- 2
""",
        )

    def test_match_score_open_penalty3(self):
        """Test 3."""
        alignments = Align.global_align(
            "GAACT", "GAT", o=-0.1, x=0  # open_gap_score, extend_gap_score
        )
        self.assertEqual(len(alignments), 1)
        self.assertAlmostEqual(alignments.score, 2.9)
        self.assertEqual(
            str(alignments[0]),
            """\
target            0 GAACT 5
                  0 ||--| 5
query             0 GA--T 3
""",
        )

    def test_match_score_open_penalty4(self):
        """Test 4."""
        alignments = Align.global_align(
            "GCT",
            "GATA",
            m=(1, -2),  # match_score, mismatch_score
            o=-0.1,  # open_gap_score
            x=0,  # extend_gap_score
        )
        self.assertEqual(len(alignments), 2)
        self.assertAlmostEqual(alignments.score, 1.7)
        self.assertEqual(
            str(alignments[0]),
            """\
target            0 G-CT- 3
                  0 |--|- 5
query             0 GA-TA 4
""",
        )
        self.assertEqual(
            str(alignments[1]),
            """\
target            0 GC-T- 3
                  0 |--|- 5
query             0 G-ATA 4
""",
        )


class TestPairwiseExtendPenalty(unittest.TestCase):
    """Alignments with gap-extend penalties."""

    def test_extend_penalty1(self):
        """Test 1."""
        alignments = Align.global_align(
            "GACT", "GT", o=-0.5, x=-0.2  # open_gap_score, extend_gap_score
        )
        self.assertEqual(len(alignments), 1)
        self.assertAlmostEqual(alignments.score, 1.3)
        self.assertEqual(
            str(alignments[0]),
            """\
target            0 GACT 4
                  0 |--| 4
query             0 G--T 2
""",
        )

    def test_extend_penalty2(self):
        """Test 2."""
        alignments = Align.global_align(
            "GACT", "GT", o=-1.5, x=-0.2  # open_gap_score, extend_gap_score
        )
        self.assertEqual(len(alignments), 1)
        self.assertAlmostEqual(alignments.score, 0.3)
        self.assertEqual(
            str(alignments[0]),
            """\
target            0 GACT 4
                  0 |--| 4
query             0 G--T 2
""",
        )


class TestPairwisePenalizeExtendWhenOpening(unittest.TestCase):
    """Alignment with ``penalize_extend_when_opening``."""

    def test_penalize_extend_when_opening(self):
        """Add gap-extend penalty to gap-opening penalty."""
        alignments = Align.global_align(
            "GACT", "GT", o=-0.2 - 1.5, x=-1.5  # open_gap_score, extend_gap_score
        )
        self.assertEqual(len(alignments), 1)
        self.assertAlmostEqual(alignments.score, -1.2)
        self.assertEqual(
            str(alignments[0]),
            """\
target            0 GACT 4
                  0 |--| 4
query             0 G--T 2
""",
        )


class TestPairwisePenalizeEndgaps(unittest.TestCase):
    """Alignments with end-gaps penalized or not."""

    def test_penalize_end_gaps(self):
        """Turn off end-gap penalties."""
        alignments = Align.global_align(
            "GACT",
            "GT",
            o=-0.8,  # open_gap_score
            x=-0.2,  # extend_gap_score
            e=0,  # end_gap_score
        )
        self.assertEqual(len(alignments), 3)
        self.assertAlmostEqual(alignments.score, 1.0)
        self.assertEqual(
            str(alignments[0]),
            """\
target            0 GACT 4
                  0 --.| 4
query             0 --GT 2
""",
        )
        self.assertEqual(
            str(alignments[1]),
            """\
target            0 GACT 4
                  0 |--| 4
query             0 G--T 2
""",
        )
        self.assertEqual(
            str(alignments[2]),
            """\
target            0 GACT 4
                  0 |.-- 4
query             0 GT-- 2
""",
        )

    def test_separate_penalize_end_gaps(self):
        """Test alignment where end-gaps are differently penalized."""
        alignments = Align.global_align(
            "AT",
            "AGG",
            m=(1.0, -0.5),  # match_score, mismatch_score
            o=-1.75,  # open_gap_score
            x=-0.25,  # extend_gap_score
            end_deletion_score=0,
        )
        self.assertAlmostEqual(alignments.score, -1.0)
        self.assertEqual(len(alignments), 1)
        self.assertEqual(
            str(alignments[0]),
            """\
target            0 A--T 2
                  0 |--- 4
query             0 AGG- 3
""",
        )


class TestPairwiseSeparateGapPenalties(unittest.TestCase):
    """Alignments with separate gap-open penalties for both sequences."""

    def test_separate_gap_penalties1(self):
        """Test 1."""
        alignments = Align.global_align(
            "GAT",
            "GTCT",
            i=(-0.3, 0),  # insertion_score
            d=(-0.8, 0),  # deletion_score
        )
        self.assertAlmostEqual(alignments.score, 1.7)
        self.assertEqual(len(alignments), 2)
        self.assertEqual(
            str(alignments[0]),
            """\
target            0 G-AT 3
                  0 |-.| 4
query             0 GTCT 4
""",
        )
        self.assertEqual(
            str(alignments[1]),
            """\
target            0 GA-T 3
                  0 |.-| 4
query             0 GTCT 4
""",
        )

    def test_separate_gap_penalties2(self):
        """Test 2."""
        alignments = Align.local_align(
            "GAT",
            "GTCT",
            i=(-0.5, 0),  # open_insertion_score, extend_insertion_score
            d=(-0.2, 0),  # open_deletion_score, extend_deletion_score
        )
        self.assertAlmostEqual(alignments.score, 1.8)
        self.assertEqual(len(alignments), 1)
        self.assertEqual(
            str(alignments[0]),
            """\
target            0 GAT 3
                  0 |-| 3
query             0 G-T 2
""",
        )


class TestPairwiseMatchDictionary(unittest.TestCase):
    """Alignments with match dictionaries."""

    match_dict = {("A", "A"): 1.5, ("A", "T"): 0.5, ("T", "T"): 1.0}
    substitution_matrix = substitution_matrices.Array(data=match_dict)

    def test_match_dictionary1(self):
        """Test 1."""
        alignments = Align.local_align(
            "ATAT",
            "ATT",
            o=-0.5,  # open_gap_score
            x=0,  # extend_gap_score
            s=self.substitution_matrix,  # substitution_matrix
        )
        self.assertAlmostEqual(alignments.score, 3.0)
        self.assertEqual(len(alignments), 2)
        self.assertEqual(
            str(alignments[0]),
            """\
target            0 ATA 3
                  0 ||. 3
query             0 ATT 3
""",
        )
        self.assertEqual(
            str(alignments[1]),
            """\
target            0 ATAT 4
                  0 ||-| 4
query             0 AT-T 3
""",
        )

    def test_match_dictionary2(self):
        """Test 2."""
        alignments = Align.local_align(
            "ATAT",
            "ATT",
            x=0,  # extend_gap_score
            s=self.substitution_matrix,  # substitution_matrix
        )
        self.assertAlmostEqual(alignments.score, 3.0)
        self.assertEqual(len(alignments), 1)
        self.assertEqual(
            str(alignments[0]),
            """\
target            0 ATA 3
                  0 ||. 3
query             0 ATT 3
""",
        )

    def test_match_dictionary3(self):
        """Test 3."""
        # Note: The PairwiseAligner distinguishes between aligning A-T and T-A.
        # The match_dict has a substitution score of 0.5 for ("A", "T") but
        # does not define a substitution score for ("T", "A"), which therefore
        # defaults to 0. On the other hand, pairwise2 uses a substitution score
        # of 0.5 both for A-T and for T-A alignments.
        alignments = Align.local_align(
            "ATT",
            "ATAT",
            o=-1,  # open_gap_score
            x=0,  # extend_gap_score
            s=self.substitution_matrix,  # substitution_matrix
        )
        self.assertAlmostEqual(alignments.score, 2.5)
        self.assertEqual(len(alignments), 2)
        self.assertEqual(
            str(alignments[0]),
            """\
target            0 AT 2
                  0 || 2
query             0 AT 2
""",
        )
        self.assertEqual(
            str(alignments[1]),
            """\
target            0 AT 2
                  0 || 2
query             2 AT 4
""",
        )


class TestPairwiseOneCharacter(unittest.TestCase):
    """Alignments where one sequence has length 1."""

    def test_align_one_char1(self):
        """Test sequence with only one match."""
        alignments = Align.local_align(
            "abcde", "c", o=-0.3, x=-0.1  # open_gap_score, extend_gap_score
        )
        self.assertAlmostEqual(alignments.score, 1.0)
        self.assertEqual(len(alignments), 1)
        self.assertEqual(
            str(alignments[0]),
            """\
target            2 c 3
                  0 | 1
query             0 c 1
""",
        )

    def test_align_one_char2(self):
        """Test sequences with two possible match positions."""
        alignments = Align.local_align(
            "abcce", "c", o=-0.3, x=-0.1  # open_gap_score, extend_gap_score
        )
        self.assertAlmostEqual(alignments.score, 1.0)
        self.assertEqual(len(alignments), 2)
        self.assertEqual(
            str(alignments[0]),
            """\
target            2 c 3
                  0 | 1
query             0 c 1
""",
        )
        self.assertEqual(
            str(alignments[1]),
            """\
target            3 c 4
                  0 | 1
query             0 c 1
""",
        )

    def test_align_one_char3(self):
        """Like test 1, but global alignment."""
        alignments = Align.global_align(
            "abcde", "c", o=-0.3, x=-0.1  # open_gap_score, extend_gap_score
        )
        self.assertAlmostEqual(alignments.score, 0.2)
        self.assertEqual(len(alignments), 1)
        self.assertEqual(
            str(alignments[0]),
            """\
target            0 abcde 5
                  0 --|-- 5
query             0 --c-- 1
""",
        )


class TestPersiteGapPenalties(unittest.TestCase):
    """Check gap penalty callbacks use correct gap opening position.

    This tests that the gap penalty callbacks are really being used
    with the correct gap opening position.
    """

    def test_gap_here_only_1(self):
        """Open a gap in second sequence only."""
        seq1 = "AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA"
        seq2 = "AABBBAAAACCCCAAAABBBAA"

        def no_gaps(x, y):
            """Very expensive to open a gap in seq1."""
            return -2000 - y

        def specific_gaps(x, y):
            """Very expensive to open a gap in seq2.

            ...unless it is in one of the allowed positions:
            """
            breaks = [0, 11, len(seq2)]
            return (-2 - y) if x in breaks else (-2000 - y)

        alignments = Align.global_align(
            seq1, seq2, i=no_gaps, d=specific_gaps  # insertion_score, deletion_score
        )
        self.assertAlmostEqual(alignments.score, 2.0)
        self.assertEqual(len(alignments), 1)
        self.assertEqual(
            str(alignments[0]),
            """\
target            0 AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA 36
                  0 --|||||||||||----------|||||||||||-- 36
query             0 --AABBBAAAACC----------CCAAAABBBAA-- 22
""",
        )

    def test_gap_here_only_2(self):
        """Force a bad alignment.

        Forces a bad alignment by having a very expensive gap penalty
        where one would normally expect a gap, and a cheap gap penalty
        in another place.
        """
        seq1 = "AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA"
        seq2 = "AABBBAAAACCCCAAAABBBAA"

        def no_gaps(x, y):
            """Very expensive to open a gap in seq1."""
            return -2000 - y

        def specific_gaps(x, y):
            """Very expensive to open a gap in seq2.

            ...unless it is in one of the allowed positions:
            """
            breaks = [0, 3, len(seq2)]
            return (-2 - y) if x in breaks else (-2000 - y)

        alignments = Align.global_align(
            seq1,
            seq2,
            m=(1, -1),  # match_score, mismatch_score
            i=no_gaps,  # insertion_score
            d=specific_gaps,  # deletion_score
        )
        self.assertAlmostEqual(alignments.score, -10.0)
        self.assertEqual(len(alignments), 2)
        self.assertEqual(
            str(alignments[0]),
            """\
target            0 AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA 36
                  0 --|||----------......|||||||||||||-- 36
query             0 --AAB----------BBAAAACCCCAAAABBBAA-- 22
""",
        )
        self.assertEqual(
            str(alignments[1]),
            """\
target            0 AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA 36
                  0 ||.------------......|||||||||||||-- 36
query             0 AAB------------BBAAAACCCCAAAABBBAA-- 22
""",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
