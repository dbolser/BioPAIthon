# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests for the validation and counting code in Bio.Align.analysis.

The numeric expectations here are worked out by hand from Nei and Gojobori
(1986) rather than taken from the output of the functions under test.
"""

import math
import unittest
import warnings

from Bio import BiopythonWarning
from Bio.Align import Alignment
from Bio.Align.analysis import _count_diff_NG86
from Bio.Align.analysis import _count_replacement
from Bio.Align.analysis import _count_site_NG86
from Bio.Align.analysis import _distance_f84
from Bio.Align.analysis import _G_test
from Bio.Align.analysis import _get_codon2codon_matrix
from Bio.Align.analysis import _get_kappa_YN00
from Bio.Align.analysis import _get_pi
from Bio.Align.analysis import _get_Q
from Bio.Align.analysis import _q
from Bio.Align.analysis import calculate_dn_ds
from Bio.Align.analysis import mktest
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

try:
    import scipy
except ImportError:
    scipy = None

CODON_TABLE = CodonTable.generic_by_id[1]


def alignment(*sequences, records=False):
    """Return an ungapped Alignment over the given sequences."""
    if records:
        return Alignment(
            [
                SeqRecord(Seq(sequence), id="seq%d" % index)
                for index, sequence in enumerate(sequences)
            ]
        )
    return Alignment([Seq(sequence) for sequence in sequences])


class CalculateDnDsValidationTests(unittest.TestCase):
    """Argument checking in calculate_dn_ds."""

    def test_cfreq_only_with_ml(self):
        """A codon frequency model is only meaningful for the ML method."""
        pair = alignment("ATGAAATTT", "ATGAAGTTT")
        for method in ("NG86", "LWL85", "YN00"):
            with self.subTest(method=method):
                with self.assertRaises(ValueError) as context:
                    calculate_dn_ds(pair, method=method, cfreq="F3x4")
                self.assertIn("ML method", str(context.exception))

    def test_unknown_cfreq(self):
        """Only F1x4, F3x4 and F61 are known codon frequency models."""
        pair = alignment("ATGAAATTT", "ATGAAGTTT")
        with self.assertRaises(ValueError) as context:
            calculate_dn_ds(pair, method="ML", cfreq="F2x4")
        self.assertIn("F1x4", str(context.exception))

    def test_unknown_method(self):
        """Only the four documented methods exist."""
        pair = alignment("ATGAAATTT", "ATGAAGTTT")
        with self.assertRaises(ValueError) as context:
            calculate_dn_ds(pair, method="NG87")
        self.assertIn("NG87", str(context.exception))

    def test_bad_character_in_target(self):
        """A non ACGT character in the first sequence is reported."""
        pair = alignment("ATGAANTTT", "ATGAAGTTT")
        with self.assertRaises(ValueError) as context:
            calculate_dn_ds(pair)
        message = str(context.exception)
        self.assertIn("AAN", message)
        self.assertIn("target", message)

    def test_bad_character_in_query(self):
        """A non ACGT character in the second sequence is reported."""
        pair = alignment("ATGAAATTT", "ATGAANTTT")
        with self.assertRaises(ValueError) as context:
            calculate_dn_ds(pair)
        message = str(context.exception)
        self.assertIn("AAN", message)
        self.assertIn("query", message)

    def test_plain_sequences_and_seqrecords_agree(self):
        """Seq and SeqRecord sequences give the same answer."""
        first = "ATGAAATTTGGGCCCAAAGGGTTTATG"
        second = "ATGAAGTTCGGACCCAAGGGGTTTATG"
        self.assertEqual(
            calculate_dn_ds(alignment(first, second)),
            calculate_dn_ds(alignment(first, second, records=True)),
        )


class NG86Tests(unittest.TestCase):
    """Nei and Gojobori (1986) site and difference counting."""

    def test_synonymous_site_counts(self):
        """Site counts follow the definition in Nei and Gojobori (1986).

        A codon has three sites; the synonymous fraction of a site is the
        fraction of the three possible changes at that position that leave
        the amino acid unchanged.  Changes to a stop codon count as
        non-synonymous.
        """
        # TTT (Phe): only TTC is synonymous, so 1 of the 9 changes.
        self.assertEqual(_count_site_NG86(["TTT"], CODON_TABLE), (1 / 3, 8 / 3))
        # GGG (Gly): all three third position changes are synonymous.
        self.assertEqual(_count_site_NG86(["GGG"], CODON_TABLE), (1.0, 2.0))
        # ATG (Met) has no synonymous changes at all.
        self.assertEqual(_count_site_NG86(["ATG"], CODON_TABLE), (0.0, 3.0))
        # TTG (Leu): CTG and TTA are synonymous, TAG is a stop codon (reached
        # by a transversion) and so counts as non-synonymous.
        self.assertEqual(_count_site_NG86(["TTG"], CODON_TABLE), (2 / 3, 7 / 3))
        # TGG (Trp) has no synonyms; two of its neighbours, TGA and TAG, are
        # stop codons reached by a transition and count as non-synonymous.
        self.assertEqual(_count_site_NG86(["TGG"], CODON_TABLE), (0.0, 3.0))
        # Counts add over the codons of a sequence
        self.assertEqual(
            _count_site_NG86(["TTT", "GGG", "ATG"], CODON_TABLE),
            (1 / 3 + 1.0, 8 / 3 + 2.0 + 3.0),
        )

    def test_transition_transversion_ratio(self):
        """k weights transitions when counting sites.

        For TTT with k = 2 the three transitions (CTT, TCT, TTC) count twice;
        only TTC is synonymous, so S = 2 and N = 2 + 2 + 6 = 10 of a total
        of 12, which normalises to 3 x 2/12 = 0.5 synonymous sites.
        """
        self.assertEqual(_count_site_NG86(["TTT"], CODON_TABLE, k=2), (0.5, 2.5))

    def test_difference_counts(self):
        """Single differences are classified as synonymous or not."""
        self.assertEqual(_count_diff_NG86("TTT", "TTT", CODON_TABLE), [0, 0])
        # TTT and TTC are both Phe
        self.assertEqual(_count_diff_NG86("TTT", "TTC", CODON_TABLE), [1, 0])
        # TTT is Phe, TTG is Leu
        self.assertEqual(_count_diff_NG86("TTT", "TTG", CODON_TABLE), [0, 1])

    def test_two_differences_average_over_both_pathways(self):
        """With two differences both mutational orders are averaged.

        TTT (Phe) to CTC (Leu) can go via CTT (Leu) or via TTC (Phe).  The
        first path is one non-synonymous change (Phe to Leu) followed by one
        synonymous change (CTT to CTC); the second is one synonymous change
        (TTT to TTC) followed by one non-synonymous one.  Each of the four
        steps carries weight 0.5, so the totals are one synonymous and one
        non-synonymous difference.
        """
        self.assertEqual(_count_diff_NG86("TTT", "CTC", CODON_TABLE), [1.0, 1.0])
        # TTA (Leu) to CTG (Leu) is synonymous by either route
        self.assertEqual(_count_diff_NG86("TTA", "CTG", CODON_TABLE), [2.0, 0.0])

    def test_identical_sequences(self):
        """Identical sequences have no substitutions at all."""
        pair = alignment("ATGAAATTTGGG", "ATGAAATTTGGG")
        self.assertEqual(calculate_dn_ds(pair, method="NG86"), (0.0, 0.0))

    def test_hand_computed_synonymous_case(self):
        """A single synonymous change, worked through by hand.

        Codons TTT GGG GGG GGG against TTC GGG GGG GGG.  Both sequences have
        1/3 + 3 = 10/3 synonymous sites and 12 - 10/3 = 26/3 non-synonymous
        sites; there is one synonymous difference and no non-synonymous
        ones.  So pS = 3/10 and, with the Jukes-Cantor correction,
        dS = -3/4 ln(1 - 4/3 pS) and dN = 0.
        """
        pair = alignment("TTTGGGGGGGGG", "TTCGGGGGGGGG")
        expected_ds = -0.75 * math.log(1 - (4 / 3) * (1 / (10 / 3)))
        dn, ds = calculate_dn_ds(pair, method="NG86")
        self.assertAlmostEqual(dn, 0.0, places=12)
        self.assertAlmostEqual(ds, expected_ds, places=12)

    def test_hand_computed_nonsynonymous_case(self):
        """A single non-synonymous change, worked through by hand.

        Codons TTT GGG GGG GGG against TTG GGG GGG GGG.  The first sequence
        has 10/3 synonymous sites and the second 2/3 + 3 = 11/3, so the
        average is 3.5 and there are 12 - 3.5 = 8.5 non-synonymous sites.
        There is one non-synonymous difference, so pN = 1/8.5 and
        dN = -3/4 ln(1 - 4/3 pN), while dS = 0.
        """
        pair = alignment("TTTGGGGGGGGG", "TTGGGGGGGGGG")
        expected_dn = -0.75 * math.log(1 - (4 / 3) * (1 / 8.5))
        dn, ds = calculate_dn_ds(pair, method="NG86")
        self.assertAlmostEqual(dn, expected_dn, places=12)
        self.assertAlmostEqual(ds, 0.0, places=12)

    def test_saturated_divergence_is_flagged(self):
        """Above 3/4 divergence the Jukes-Cantor correction is undefined.

        The implementation reports -1 instead of raising, both for dN and
        for dS.
        """
        # every codon differs at every position and nothing is synonymous
        pair = alignment("TTTTTTTTT", "GGGGGGGGG")
        dn, ds = calculate_dn_ds(pair, method="NG86")
        self.assertEqual(dn, -1)
        # a single synonymous difference in a single codon saturates pS
        pair = alignment("TTT", "TTC")
        dn, ds = calculate_dn_ds(pair, method="NG86")
        self.assertEqual(ds, -1)

    def test_dn_ds_rises_with_divergence(self):
        """More differences means larger distances."""
        reference = "ATGAAATTTGGGCCCAAAGGGTTTATGCCCAAAGGG"
        one = "ATGAAGTTTGGGCCCAAAGGGTTTATGCCCAAAGGG"
        two = "ATGAAGTTCGGACCCAAGGGGTTTATGCCCAAAGGG"
        distances = [
            calculate_dn_ds(alignment(reference, other), method="NG86")[1]
            for other in (reference, one, two)
        ]
        self.assertEqual(distances, sorted(distances))


class DistanceF84Tests(unittest.TestCase):
    """The F84 distance, with PAML's fallbacks to K80 and to JC69."""

    equal = {"T": 0.25, "C": 0.25, "A": 0.25, "G": 0.25}

    def k80(self, P, Q):
        """Return kappa and the distance under K80 (Kimura 1980)."""
        a = -math.log(1 - 2 * P - Q)
        b = -math.log(1 - 2 * Q)
        return (a / 2 - b / 4) / (b / 4), a / 2 + b / 4

    def test_equal_base_frequencies_reduce_to_k80(self):
        """With equal base frequencies F84 is K80."""
        kappa, t = _distance_f84(100, 0.1, 0.05, self.equal)
        expected_kappa, expected_t = self.k80(0.1, 0.05)
        self.assertAlmostEqual(kappa, expected_kappa, places=12)
        self.assertAlmostEqual(t, expected_t, places=12)

    def test_k80_when_f84_fails(self):
        """With no pyrimidines F84 is undefined, so K80 is used."""
        pi = {"T": 0, "C": 0, "A": 0.5, "G": 0.5}
        self.assertEqual(_distance_f84(100, 0.1, 0.05, pi), self.k80(0.1, 0.05))

    def test_jc69_when_there_are_no_transversions(self):
        """Without transversions neither F84 nor K80 gives a kappa."""
        kappa, t = _distance_f84(100, 0.1, 0, self.equal)
        self.assertEqual(kappa, -1)
        self.assertAlmostEqual(t, -0.75 * math.log(1 - 4 / 3 * 0.1), places=12)
        # yn00.c counts a proportion below 1e-10 as no transversions
        kappa, t = _distance_f84(100, 0.1, 1e-11, self.equal)
        self.assertEqual(kappa, -1)
        self.assertAlmostEqual(t, -0.75 * math.log(1 - 4 / 3 * 0.1), places=10)

    def test_maximum_kappa(self):
        """Kappa is at most 999."""
        kappa, t = _distance_f84(100, 0.1, 1e-8, self.equal)
        self.assertEqual(kappa, 999)

    def test_no_differences(self):
        """No differences means a distance of zero, and not of -0.0."""
        kappa, t = _distance_f84(100, 0, 0, self.equal)
        self.assertEqual(kappa, -1)
        self.assertEqual(t, 0)
        self.assertEqual(math.copysign(1, t), 1)

    def test_saturated(self):
        """Beyond what JC69 can correct, the distance is 3/4 log(n), at most 99."""
        # JC69 needs P + Q < 3/4
        kappa, t = _distance_f84(10, 0.8, 0, self.equal)
        self.assertEqual(kappa, -1)
        self.assertAlmostEqual(t, 0.75 * math.log(10), places=12)
        # yn00.c would give a distance of zero or less for one site or fewer
        self.assertEqual(_distance_f84(0.5, 0.8, 0, self.equal), (-1, 99))
        # 3/4 log(n) is more than 99 for n above about 2e57
        self.assertEqual(_distance_f84(1e60, 0.8, 0, self.equal), (-1, 99))
        # more differences than sites
        self.assertEqual(_distance_f84(0.5, 0.6, 0.6, self.equal), (1, 99))


@unittest.skipIf(scipy is None, "the YN00 method needs scipy")
class YN00Tests(unittest.TestCase):
    """The Yang and Nielsen (2000) method.

    Unless a test says otherwise, the expected values are those of the yn00
    program of PAML 4.10.10, run with weighting = 1, to the four decimal
    places it prints. yn00 gives the same values with the two sequences in
    either order.
    """

    def yn00(self, first, second, codon_table=CODON_TABLE):
        """Return dN and dS by YN00 with the sequences both ways round.

        Fails on any warning, and if the order of the sequences matters.
        """
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            dn, ds = calculate_dn_ds(
                alignment(first, second), method="YN00", codon_table=codon_table
            )
            reverse = calculate_dn_ds(
                alignment(second, first), method="YN00", codon_table=codon_table
            )
        self.assertAlmostEqual(reverse[0], dn, places=12)
        self.assertAlmostEqual(reverse[1], ds, places=12)
        return dn, ds

    def test_tutorial_alignment(self):
        """The E. coli and Shigella codon alignment of the Tutorial.

        This returned None. The base frequencies at the 0-fold sites were not
        all normalised, giving a kappa of -35, and every path between TAT and
        ATC went through a codon absent from both sequences, which had a
        frequency of zero, so its probability was 0/0.
        """
        dn, ds = self.yn00(
            "TCCGACGTTCTTAGGGGGTATCTTGACTACGACGCCAAGAAAGAGCGT",
            "TCCGACGTTCTTAGGGGTATCTTTGACTACGACGCCAAGAAAGAGCGT",
        )
        self.assertAlmostEqual(dn, 0.0785, places=4)
        self.assertAlmostEqual(ds, 0.3384, places=4)

    def test_issue_1611(self):
        """The pair from biopython/biopython#1611 raised ZeroDivisionError."""
        dn, ds = self.yn00("ATGGGTACCGTGTTGTAT", "ATGGCTACCGTGTTAGTA")
        self.assertAlmostEqual(dn, 0.3123, places=4)
        self.assertAlmostEqual(ds, 0.5449, places=4)

    def test_order_of_the_sequences(self):
        """Paths with two or three differences are counted step by step.

        The step to the second sequence's codon was counted as a second copy
        of the first step, and on three-difference paths it was classed at
        the position of the middle step, so the result depended on which
        sequence came first.
        """
        dn, ds = self.yn00("GATGATCGTACGGCC", "GTTTATGGTGTTACT")
        self.assertAlmostEqual(dn, 0.8158, places=4)
        self.assertAlmostEqual(ds, 0.6176, places=4)

    def test_identical_sequences(self):
        """Identical sequences raised "dictionary changed size"."""
        dn, ds = self.yn00("ATGAAATTTGGG", "ATGAAATTTGGG")
        self.assertEqual((dn, ds), (0, 0))
        self.assertEqual((math.copysign(1, dn), math.copysign(1, ds)), (1, 1))

    def test_synonymous_differences_only(self):
        """Only synonymous differences: dN is zero."""
        dn, ds = self.yn00("ATGAAATTTGGGCTG", "ATGAAGTTCGGACTG")
        self.assertEqual(dn, 0)
        self.assertAlmostEqual(ds, 2.3933, places=4)

    def test_nonsynonymous_differences_only(self):
        """Only nonsynonymous differences: dS is zero.

        dN was 0.3236 when the base frequencies at the nonsynonymous sites
        counted the base of each neighbouring codon, not the codon's own.
        """
        dn, ds = self.yn00("ATGAAATTTGGG", "ATGAACTTAGGG")
        self.assertAlmostEqual(dn, 0.3255, places=4)
        self.assertEqual(ds, 0)

    def test_every_path_through_a_stop_codon(self):
        """Where every path passes through a stop codon, count as yn00 does.

        In the vertebrate mitochondrial code AGA, AGG and TAG are stop
        codons, so every path from TGG to AAA passes through one. yn00
        counts the three differences as nonsynonymous, half a transition
        and two and a half transversions.
        """
        table = CodonTable.generic_by_id[2]
        dn, ds = self.yn00("ATGTGGCTGAAACCC", "ATGAAACTGAAACCC", table)
        self.assertAlmostEqual(dn, 0.5398, places=4)
        self.assertEqual(ds, 0)

    def test_kappa_where_k80_fails(self):
        """A class of sites whose K80 distance is undefined gives no kappa.

        At the three fourfold degenerate third positions here, 1 - 2P - Q is
        exactly 0, so only the non-degenerate sites give an estimate. P and
        Q have to be computed as yn00 computes them: as differences over
        sites, 1 - 2P - Q rounds to just above 0, and kappa comes out at 13.3.
        """
        codons1 = "TCT AGC AAG CTT GGG CAA AAG GTT TCA CAC".split()
        codons2 = "TCC GGT AAG GTT ACG CAG GGC GTT TCT CAC".split()
        kappa = _get_kappa_YN00(codons1, codons2, CODON_TABLE)
        self.assertAlmostEqual(kappa, 5.1727, places=4)

    def test_one_sequence_without_synonymous_sites(self):
        """The base frequencies at synonymous sites come from the other.

        Methionine and tryptophan have one codon each, so ATG and TGG have no
        synonymous sites. yn00 averages the base frequencies at the
        synonymous sites of the two sequences, and as those of the second
        are 0/0, it gives nan. These values are this implementation's own.
        """
        dn, ds = self.yn00("TGCCTGCGTATAAGG", "TGGATGTGGATGTGG")
        self.assertAlmostEqual(dn, 0.8334, places=4)
        self.assertAlmostEqual(ds, 0.0634, places=4)

    def test_no_synonymous_or_no_nonsynonymous_sites(self):
        """Without both kinds of sites, dN and dS are nan, with a warning.

        These raised ZeroDivisionError. yn00 gives nan for omega and for at
        least one of dN and dS.
        """
        for first, second in (
            ("ATGTGGATG", "ATGTGGTGG"),  # no synonymous sites
            ("CAGCAGCAGCAG", "CAGCAACAGCAG"),  # no nonsynonymous sites
            ("AAAAAAAAA", "AAAAAAAAA"),  # no sites
        ):
            with self.subTest(first=first, second=second):
                pair = alignment(first, second)
                with self.assertWarnsRegex(
                    BiopythonWarning, "needs synonymous and nonsynonymous sites"
                ):
                    dn, ds = calculate_dn_ds(pair, method="YN00")
                self.assertTrue(math.isnan(dn))
                self.assertTrue(math.isnan(ds))

    def test_saturated_and_not_converging(self):
        """A saturated pair gives yn00's maximum dS, and a warning.

        The estimates of dN alternate between 1.6102 and 2.8260, as they do
        in yn00, and the last ones are returned, as yn00 returns them, with a
        warning.
        """
        pair = alignment("TTTTTTTTT", "CCACCACCA")
        with self.assertWarnsRegex(BiopythonWarning, "did not converge"):
            dn, ds = calculate_dn_ds(pair, method="YN00")
        self.assertAlmostEqual(dn, 2.8260, places=4)
        self.assertEqual(ds, 99)


class GetPiTests(unittest.TestCase):
    """Codon frequency models used by the ML method."""

    codons1 = ["ATG", "AAA", "TTT"]
    codons2 = ["ATG", "AAG", "TTC"]

    def test_f3x4(self):
        """F3x4 multiplies the base frequencies of the three positions."""
        pi = _get_pi(self.codons1, self.codons2, "F3x4", CODON_TABLE)
        codons = self.codons1 + self.codons2
        frequencies = []
        for position in range(3):
            bases = [codon[position] for codon in codons]
            frequencies.append(
                {base: bases.count(base) / len(bases) for base in "ACGT"}
            )
        for codon in ("ATG", "AAA", "TTT", "TAA"):
            with self.subTest(codon=codon):
                expected = (
                    frequencies[0][codon[0]]
                    * frequencies[1][codon[1]]
                    * frequencies[2][codon[2]]
                )
                self.assertAlmostEqual(pi[codon], expected, places=12)
        # every sense codon and every stop codon is present
        self.assertEqual(len(pi), 64)

    def test_f1x4(self):
        """F1x4 uses one set of base frequencies for all three positions."""
        pi = _get_pi(self.codons1, self.codons2, "F1x4", CODON_TABLE)
        codons = self.codons1 + self.codons2
        bases = "".join(codons)
        frequency = {base: bases.count(base) / len(bases) for base in "ACGT"}
        for codon in ("ATG", "AAA", "TTT"):
            expected = frequency[codon[0]] * frequency[codon[1]] * frequency[codon[2]]
            self.assertAlmostEqual(pi[codon], expected, places=12)

    def test_f1x4_base_absent_from_both_sequences(self):
        """F1x4 gives a base that does not occur a frequency of zero."""
        # There is no C. This used to raise KeyError, although F3x4 already
        # gave such codons zero, as did the copy of this code in codonalign.
        pi = _get_pi(["ATG", "AAA"], ["ATG", "AAG"], "F1x4", CODON_TABLE)
        self.assertEqual(len(pi), 64)
        self.assertEqual(pi["CCC"], 0)
        self.assertEqual(pi["ACG"], 0)
        # A, T and G occur 7, 2 and 3 times in 12 bases
        self.assertAlmostEqual(pi["ATG"], 7 / 12 * 2 / 12 * 3 / 12, places=12)

    def test_f61(self):
        """F61 counts whole codons, with a pseudo count of 0.1."""
        pi = _get_pi(self.codons1, self.codons2, "F61", CODON_TABLE)
        total = 64 * 0.1 + len(self.codons1) + len(self.codons2)
        self.assertAlmostEqual(pi["ATG"], (0.1 + 2) / total, places=12)
        self.assertAlmostEqual(pi["AAA"], (0.1 + 1) / total, places=12)
        self.assertAlmostEqual(pi["GGG"], 0.1 / total, places=12)


class QMatrixTests(unittest.TestCase):
    """The codon substitution rate matrix used by the ML method."""

    pi = dict.fromkeys(
        [
            first + second + third
            for first in "ACGT"
            for second in "ACGT"
            for third in "ACGT"
        ],
        1 / 64,
    )

    def rate(self, codon1, codon2, k=2.0, w=0.5):
        """Return the off diagonal rate for one codon pair."""
        return _q(codon1, codon2, self.pi, k, w, CODON_TABLE)

    def test_diagonal_is_zero(self):
        """A codon does not substitute for itself."""
        self.assertEqual(self.rate("ATG", "ATG"), 0)

    def test_stop_codons_are_excluded(self):
        """No substitution involves a stop codon."""
        self.assertEqual(self.rate("TAA", "TAC"), 0)
        self.assertEqual(self.rate("TAC", "TAA"), 0)

    def test_codons_missing_from_pi_are_excluded(self):
        """A codon with no frequency contributes no rate."""
        pi = dict(self.pi)
        del pi["TTC"]
        self.assertEqual(_q("TTT", "TTC", pi, 2.0, 0.5, CODON_TABLE), 0)
        self.assertEqual(_q("TTC", "TTT", pi, 2.0, 0.5, CODON_TABLE), 0)

    def test_multiple_differences_are_excluded(self):
        """Only single nucleotide substitutions have a rate."""
        self.assertEqual(self.rate("TTT", "CTC"), 0)
        self.assertEqual(self.rate("TTT", "GGG"), 0)

    def test_rate_categories(self):
        """The four rate categories differ by factors of k and w.

        Following Goldman and Yang (1994), a synonymous transversion has
        rate pi_j, a synonymous transition k pi_j, a non-synonymous
        transversion w pi_j and a non-synonymous transition w k pi_j.
        """
        k, w, frequency = 2.0, 0.5, 1 / 64
        # TTT -> TTC: Phe to Phe, pyrimidine transition
        self.assertAlmostEqual(self.rate("TTT", "TTC"), k * frequency)
        # AAA -> AAG: Lys to Lys, purine transition
        self.assertAlmostEqual(self.rate("AAA", "AAG"), k * frequency)
        # CTT -> CTA: Leu to Leu, pyrimidine to purine transversion
        self.assertAlmostEqual(self.rate("CTT", "CTA"), frequency)
        # AAA -> GAA: Lys to Glu, purine transition
        self.assertAlmostEqual(self.rate("AAA", "GAA"), w * k * frequency)
        # TTT -> TTA: Phe to Leu, transversion
        self.assertAlmostEqual(self.rate("TTT", "TTA"), w * frequency)

    def test_get_Q_rows_sum_to_zero(self):
        """Each row of the rate matrix sums to zero."""
        codons = ["TTT", "TTC", "TTA", "TTG", "CTT", "ATG"]
        matrix = _get_Q(self.pi, 2.0, 0.5, codons, CODON_TABLE)
        for row in matrix:
            self.assertAlmostEqual(sum(row), 0.0, places=12)

    def test_get_Q_tolerates_codons_missing_from_pi(self):
        """A codon absent from pi does not break the normalisation."""
        pi = dict(self.pi)
        del pi["ATG"]
        codons = ["TTT", "TTC", "TTA", "TTG", "CTT", "ATG"]
        matrix = _get_Q(pi, 2.0, 0.5, codons, CODON_TABLE)
        for row in matrix:
            self.assertAlmostEqual(sum(row), 0.0, places=12)


class CountReplacementTests(unittest.TestCase):
    """Counting the substitutions needed to explain a set of codons."""

    @classmethod
    def setUpClass(cls):
        """Build the codon to codon distance matrices once."""
        cls.G, cls.nonsyn_G = _get_codon2codon_matrix(codon_table=CODON_TABLE)

    def test_single_codon_needs_no_change(self):
        """One codon needs no substitutions of either kind."""
        self.assertEqual(_count_replacement({"ATG"}, self.G), (0, 0))
        self.assertEqual(_count_replacement({"ATG"}, self.nonsyn_G), (0, 0))

    def test_two_codons(self):
        """Two codons need the substitutions that separate them."""
        # TTT and TTC differ at one position and encode the same amino acid
        self.assertEqual(_count_replacement({"TTT", "TTC"}, self.G), 1)
        self.assertEqual(_count_replacement({"TTT", "TTC"}, self.nonsyn_G), 0)
        # TTT and TTA differ at one position and change the amino acid
        self.assertEqual(_count_replacement({"TTT", "TTA"}, self.nonsyn_G), 1)

    def test_three_codons_use_a_spanning_tree(self):
        """More than two codons should be joined by a minimum spanning tree.

        TTT, TTC and TTA are mutually one substitution apart, so the minimum
        spanning tree over them costs 2.
        """
        codons = {"TTT", "TTC", "TTA"}
        self.assertEqual(_count_replacement(codons, self.G), 2)


class MkTestTests(unittest.TestCase):
    """The McDonald-Kreitman test."""

    # Two species, two sequences each, four codons per sequence:
    #  1. TTT/TTT vs TTC/TTC  - fixed between species, synonymous
    #  2. TTT/TTT vs TTA/TTA  - fixed between species, non-synonymous
    #  3. TTT/TTC vs TTT/TTT  - polymorphic within species A, synonymous
    #  4. TTT/TTA vs TTT/TTT  - polymorphic within species A, non-synonymous
    # so the contingency table is [1, 1, 1, 1].
    sequences = [
        "TTTTTTTTTTTT",
        "TTTTTTTTCTTA",
        "TTCTTATTTTTT",
        "TTCTTATTTTTT",
    ]
    species = ["A", "A", "B", "B"]

    @staticmethod
    def g_test(counts):
        """Return the p value of a 2x2 G test, from the published formula.

        G = 2 sum(O ln(O/E)) is compared against a chi-square distribution
        with one degree of freedom, whose upper tail probability is
        erfc(sqrt(G/2)).
        """
        total = sum(counts)
        rows = (counts[0] + counts[1], counts[2] + counts[3])
        columns = (counts[0] + counts[2], counts[1] + counts[3])
        expected = [
            rows[0] * columns[0] / total,
            rows[0] * columns[1] / total,
            rows[1] * columns[0] / total,
            rows[1] * columns[1] / total,
        ]
        g = 2 * sum(
            observed * math.log(observed / value)
            for observed, value in zip(counts, expected)
            if observed
        )
        return math.erfc(math.sqrt(g / 2))

    def test_balanced_table_is_not_significant(self):
        """A table of four equal counts gives a p value of exactly one."""
        self.assertAlmostEqual(
            mktest(alignment(*self.sequences), self.species),
            self.g_test([1, 1, 1, 1]),
            places=12,
        )
        self.assertAlmostEqual(
            mktest(alignment(*self.sequences), self.species), 1.0, places=12
        )

    def test_unbalanced_table(self):
        """An extra fixed synonymous site gives the table [2, 1, 1, 1]."""
        sequences = [
            sequence + extra
            for sequence, extra in zip(self.sequences, ("TTT", "TTT", "TTC", "TTC"))
        ]
        self.assertAlmostEqual(
            mktest(alignment(*sequences), self.species),
            self.g_test([2, 1, 1, 1]),
            places=12,
        )

    def test_plain_sequences_and_seqrecords_agree(self):
        """Seq and SeqRecord sequences give the same answer."""
        self.assertEqual(
            mktest(alignment(*self.sequences), self.species),
            mktest(alignment(*self.sequences, records=True), self.species),
        )

    def test_invariant_sites_are_skipped(self):
        """Codons that are identical everywhere contribute nothing."""
        padded = [sequence + "GGGGGG" for sequence in self.sequences]
        self.assertAlmostEqual(
            mktest(alignment(*padded), self.species),
            mktest(alignment(*self.sequences), self.species),
            places=12,
        )

    def test_identical_sequences(self):
        """An alignment with nothing to count is not evidence of anything.

        Every codon is invariant, so the contingency table is [0, 0, 0, 0]
        and the p value is one.
        """
        sequences = ["TTTTTTTTTTTT"] * 4
        self.assertEqual(mktest(alignment(*sequences), self.species), 1.0)

    def test_no_polymorphism(self):
        """A table with an empty row is degenerate and gives a p value of one.

        Both differences are fixed between the species, so the table is
        [1, 1, 0, 0]: the polymorphic row is empty, the observed counts equal
        the expected counts exactly, and G is zero.
        """
        sequences = ["TTTTTT", "TTTTTT", "TTCTTA", "TTCTTA"]
        self.assertEqual(mktest(alignment(*sequences), self.species), 1.0)


class GTestTests(unittest.TestCase):
    """Empty cells in the 2x2 contingency table of the G test."""

    def test_published_example(self):
        """The worked example from the docstring of _G_test."""
        # tot = 68, row totals 24 and 44, column totals 59 and 9, so the
        # expected counts are 24*59/68, 24*9/68, 44*59/68 and 44*9/68.
        counts = [17, 7, 42, 2]
        expected = [24 * 59 / 68, 24 * 9 / 68, 44 * 59 / 68, 44 * 9 / 68]
        g = 2 * sum(
            observed * math.log(observed / value)
            for observed, value in zip(counts, expected)
        )
        self.assertAlmostEqual(_G_test(counts), math.erfc(math.sqrt(g / 2)), places=12)

    def test_one_empty_cell(self):
        """A zero observed count contributes nothing to G.

        Following ``lim x->0 of x ln(x) = 0``, the empty cell drops out of
        the sum rather than raising.
        """
        # tot = 6, row totals 2 and 4, column totals 3 and 3, so the expected
        # counts are exactly [1, 1, 2, 2] and only three cells contribute.
        counts = [2, 0, 1, 3]
        g = 2 * (2 * math.log(2 / 1) + 1 * math.log(1 / 2) + 3 * math.log(3 / 2))
        self.assertAlmostEqual(_G_test(counts), math.erfc(math.sqrt(g / 2)), places=12)

    def test_empty_row(self):
        """With no polymorphism the observed counts are the expected counts."""
        # tot = 8, row totals 8 and 0, column totals 3 and 5, so the expected
        # counts are [3, 5, 0, 0], which is the table itself, and G is zero.
        self.assertEqual(_G_test([3, 5, 0, 0]), 1.0)
        # the same holds for an empty row of fixed differences
        self.assertEqual(_G_test([0, 0, 3, 5]), 1.0)

    def test_empty_column(self):
        """With no non-synonymous change the table is degenerate as well."""
        # tot = 8, row totals 3 and 5, column totals 8 and 0, so the expected
        # counts are [3, 0, 5, 0], which is the table itself, and G is zero.
        self.assertEqual(_G_test([3, 0, 5, 0]), 1.0)
        self.assertEqual(_G_test([0, 3, 0, 5]), 1.0)

    def test_empty_table(self):
        """An all-zero table has no expected counts and no information."""
        self.assertEqual(_G_test([0, 0, 0, 0]), 1.0)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
