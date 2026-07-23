# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests for the Bio.SeqUtils.MeltingTemp module.

The expected values in this file are, wherever a number is checked, derived
from the published formulae rather than from the output of the functions
under test.  The reference implementations at the top of the file are
deliberately written from the papers cited in ``Bio.SeqUtils.MeltingTemp``
so that a wrong answer from the module is not silently blessed.
"""

import math
import unittest
import warnings

from Bio import BiopythonWarning
from Bio.Seq import reverse_complement
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

# Universal gas constant in cal/(K mol), as used by the nearest neighbour
# model of SantaLucia and co-workers.
GAS_CONSTANT = 1.987


def reference_salt_correction(Na=0, K=0, Tris=0, Mg=0, dNTPs=0, method=1, seq=None):
    """Return the salt correction term, computed from the published formulae.

    The sodium equivalent conversion is von Ahsen et al. (2001), Clin Chem
    47: 1956-1961:

        [Na_eq] = [Na+] + [K+] + [Tris]/2 + 120 * sqrt([Mg2+] - [dNTPs])

    and the seven corrections are the ones listed in the module docstring.
    All concentrations are given in millimolar and converted to molar before
    being used in the logarithms.
    """
    if not method:
        return 0
    monovalent = Na + K + Tris / 2.0
    if sum((K, Mg, Tris, dNTPs)) > 0 and method != 7 and dNTPs < Mg:
        monovalent += 120 * math.sqrt(Mg - dNTPs)
    mon = monovalent * 1e-3
    if method == 1:
        return 16.6 * math.log10(mon)
    if method == 2:
        return 16.6 * math.log10(mon / (1.0 + 0.7 * mon))
    if method == 3:
        return 12.5 * math.log10(mon)
    if method == 4:
        return 11.7 * math.log10(mon)
    gc_fraction = (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0
    if method == 5:
        return 0.368 * (len(seq) - 1) * math.log(mon)
    if method == 6:
        return (4.29 * gc_fraction - 3.95) * 1e-5 * math.log(mon) + (
            9.40e-6 * math.log(mon) ** 2
        )
    # Owczarzy et al. (2008), Biochemistry 47: 5336-5353.
    a, b, c, d, e, f, g = 3.92, -0.911, 6.26, 1.42, -48.2, 52.5, 8.31
    mg = Mg * 1e-3
    if dNTPs > 0:
        dntps = dNTPs * 1e-3
        ka = 3e4  # association constant for Mg:dNTP
        mg = (
            -(ka * dntps - ka * mg + 1.0)
            + math.sqrt((ka * dntps - ka * mg + 1.0) ** 2 + 4.0 * ka * mg)
        ) / (2.0 * ka)
    if monovalent > 0:
        ratio = math.sqrt(mg) / mon
        if ratio < 0.22:
            # monovalent ions dominate; use the Owczarzy 2004 correction
            return (4.29 * gc_fraction - 3.95) * 1e-5 * math.log(mon) + (
                9.40e-6 * math.log(mon) ** 2
            )
        if ratio < 6.0:
            a = 3.92 * (0.843 - 0.352 * math.sqrt(mon) * math.log(mon))
            d = 1.42 * (1.279 - 4.03e-3 * math.log(mon) - 8.03e-3 * math.log(mon) ** 2)
            g = 8.31 * (0.486 - 0.258 * math.log(mon) + 5.25e-3 * math.log(mon) ** 3)
    return (
        a
        + b * math.log(mg)
        + gc_fraction * (c + d * math.log(mg))
        + (1 / (2.0 * (len(seq) - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
    ) * 1e-5


def reference_initiation(seq, nn_table):
    """Return the initiation (dH, dS) of a duplex, from the NN model.

    The nearest neighbour model adds a duplex initiation term, a term that
    depends on whether the duplex contains any G or C, a penalty when the 5'
    end is a T (or, symmetrically, the 3' end an A), and a per-end term that
    differs for A/T and G/C terminal base pairs.
    """
    delta_h, delta_s = nn_table["init"]
    if set(seq) <= set("AT"):
        extra_h, extra_s = nn_table["init_allA/T"]
    else:
        extra_h, extra_s = nn_table["init_oneG/C"]
    delta_h += extra_h
    delta_s += extra_s
    for applies in (seq.startswith("T"), seq.endswith("A")):
        if applies:
            delta_h += nn_table["init_5T/A"][0]
            delta_s += nn_table["init_5T/A"][1]
    for end in (seq[0], seq[-1]):
        key = "init_A/T" if end in "AT" else "init_G/C"
        delta_h += nn_table[key][0]
        delta_s += nn_table[key][1]
    return delta_h, delta_s


def reference_stacking(seq, c_seq, nn_table, imm_table):
    """Return the (dH, dS) summed over every nearest neighbour of a duplex.

    ``seq`` and ``c_seq`` must have the same length and be aligned, with
    ``c_seq`` running 3'->5'.  A dimer may be listed in the tables either as
    written or read from the other strand, which is the same duplex reversed.
    """
    delta_h = delta_s = 0.0
    for index in range(len(seq) - 1):
        key = seq[index : index + 2] + "/" + c_seq[index : index + 2]
        for table in (imm_table, nn_table):
            if table is None:
                continue
            if key in table:
                value = table[key]
                break
            if key[::-1] in table:
                value = table[key[::-1]]
                break
        else:
            raise AssertionError(f"no thermodynamic data for {key}")
        delta_h += value[0]
        delta_s += value[1]
    return delta_h, delta_s


def reference_tm_nn(
    seq,
    c_seq=None,
    nn_table=None,
    imm_table=None,
    init_seq=None,
    extra=(0.0, 0.0),
    dnac1=25,
    dnac2=25,
    selfcomp=False,
):
    """Return the nearest neighbour Tm of an aligned duplex.

    Tm = 1000 * dH / (dS + R ln k) - 273.15, with k the effective strand
    concentration: (dnac1 - dnac2/2) for a non-self-complementary duplex and
    dnac1 for a self-complementary one, which additionally carries the
    symmetry term of the table.  ``extra`` allows the caller to add the
    dangling end or terminal mismatch contributions of a duplex whose
    stacking region is shorter than the full oligonucleotide.
    """
    if nn_table is None:
        nn_table = mt.DNA_NN3
    if c_seq is None:
        c_seq = str(Seq(seq).complement())
    if init_seq is None:
        init_seq = seq
    delta_h, delta_s = reference_initiation(init_seq, nn_table)
    stack_h, stack_s = reference_stacking(seq, c_seq, nn_table, imm_table)
    delta_h += stack_h + extra[0]
    delta_s += stack_s + extra[1]
    if selfcomp:
        k = dnac1 * 1e-9
        delta_h += nn_table["sym"][0]
        delta_s += nn_table["sym"][1]
    else:
        k = (dnac1 - dnac2 / 2.0) * 1e-9
    return (1000 * delta_h) / (delta_s + GAS_CONSTANT * math.log(k)) - 273.15


class SaltCorrectionTests(unittest.TestCase):
    """Tests for Bio.SeqUtils.MeltingTemp.salt_correction."""

    conditions = [
        {"Na": 50},
        {"Na": 100},
        {"Na": 1000},
        {"Na": 50, "K": 50},
        {"Na": 100, "Tris": 20},
        {"Na": 100, "Tris": 20, "Mg": 1.5},
        {"Na": 100, "Tris": 20, "Mg": 1.5, "dNTPs": 0.6},
        {"Na": 50, "Mg": 3.0, "dNTPs": 0.8},
        {"K": 50, "Tris": 10},
        # dNTPs >= Mg, so free magnesium is ignored
        {"Na": 50, "Mg": 1.0, "dNTPs": 1.0},
        {"Na": 50, "Mg": 1.0, "dNTPs": 2.0},
    ]
    seq = "CGTTCCAAAGATGTGGGCATGAGCTTAC"

    def test_matches_published_formulae(self):
        """Every method agrees with an independent implementation."""
        for method in range(1, 8):
            for condition in self.conditions:
                with self.subTest(method=method, **condition):
                    expected = reference_salt_correction(
                        method=method, seq=self.seq, **condition
                    )
                    self.assertAlmostEqual(
                        mt.salt_correction(method=method, seq=self.seq, **condition),
                        expected,
                        places=10,
                    )

    def test_method_zero_is_no_correction(self):
        """Method 0 (or None) returns no correction at all."""
        self.assertEqual(mt.salt_correction(Na=50, method=0), 0)
        self.assertEqual(mt.salt_correction(Na=50, method=None), 0)

    def test_relative_scaling_of_the_log_methods(self):
        """Methods 1, 3 and 4 differ only in their leading constant."""
        method1 = mt.salt_correction(Na=50, method=1)
        self.assertAlmostEqual(
            mt.salt_correction(Na=50, method=3), method1 * 12.5 / 16.6, places=10
        )
        self.assertAlmostEqual(
            mt.salt_correction(Na=50, method=4), method1 * 11.7 / 16.6, places=10
        )

    def test_correction_increases_with_salt(self):
        """More counter ions means a smaller (less negative) correction."""
        for method in (1, 2, 3, 4):
            with self.subTest(method=method):
                values = [
                    mt.salt_correction(Na=concentration, method=method)
                    for concentration in (10, 50, 100, 500, 900)
                ]
                self.assertEqual(values, sorted(values))
                self.assertLess(values[-1], 0)  # below 1 M the term is negative
                if method != 2:
                    # at exactly 1 M sodium the plain log10 corrections vanish
                    self.assertAlmostEqual(
                        mt.salt_correction(Na=1000, method=method), 0, places=10
                    )

    def test_method_seven_without_magnesium_matches_method_six(self):
        """Without Mg2+ the Owczarzy 2008 tree falls back on Owczarzy 2004.

        With [Mg2+] = 0 the ratio sqrt([Mg2+])/[Mon+] is zero, which is below
        the 0.22 threshold of the decision tree, so method 7 must return
        exactly what method 6 returns.
        """
        for sodium in (10, 50, 200, 1000):
            with self.subTest(Na=sodium):
                self.assertAlmostEqual(
                    mt.salt_correction(Na=sodium, method=7, seq=self.seq),
                    mt.salt_correction(Na=sodium, method=6, seq=self.seq),
                    places=12,
                )

    def test_sodium_equivalent_conversion(self):
        """Divalent and other ions are folded into a sodium equivalent.

        von Ahsen et al. (2001): [Na_eq] = [Na+] + [K+] + [Tris]/2 +
        120 * sqrt([Mg2+] - [dNTPs]).
        """
        equivalent = 100 + 25 + 20 / 2 + 120 * math.sqrt(1.5 - 0.5)
        self.assertAlmostEqual(
            mt.salt_correction(Na=100, K=25, Tris=20, Mg=1.5, dNTPs=0.5, method=1),
            16.6 * math.log10(equivalent * 1e-3),
            places=10,
        )
        # If [dNTPs] >= [Mg2+] the magnesium term is dropped entirely
        self.assertAlmostEqual(
            mt.salt_correction(Na=100, K=25, Tris=20, Mg=1.5, dNTPs=1.5, method=1),
            16.6 * math.log10((100 + 25 + 10) * 1e-3),
            places=10,
        )

    def test_sequence_required_for_methods_five_six_and_seven(self):
        """Methods that need the sequence complain when it is missing."""
        for method in (5, 6, 7):
            with self.subTest(method=method):
                with self.assertRaises(ValueError) as context:
                    mt.salt_correction(Na=50, method=method)
                self.assertIn("sequence is missing", str(context.exception))

    def test_zero_ion_concentration_rejected(self):
        """A total ion concentration of zero has no logarithm."""
        for method in range(1, 7):
            with self.subTest(method=method):
                with self.assertRaises(ValueError) as context:
                    mt.salt_correction(method=method, seq=self.seq)
                self.assertIn("zero", str(context.exception))

    def test_unknown_method_rejected(self):
        """Only methods 1 to 7 exist."""
        with self.assertRaises(ValueError) as context:
            mt.salt_correction(Na=50, method=8, seq=self.seq)
        self.assertIn("1-7", str(context.exception))


class ChemCorrectionTests(unittest.TestCase):
    """Tests for Bio.SeqUtils.MeltingTemp.chem_correction."""

    def test_no_additives(self):
        """Without additives the melting temperature is returned unchanged."""
        self.assertEqual(mt.chem_correction(70), 70)

    def test_dmso(self):
        """DMSO lowers Tm by DMSOfactor degrees per percent."""
        self.assertAlmostEqual(mt.chem_correction(70, DMSO=3), 70 - 0.75 * 3)
        self.assertAlmostEqual(
            mt.chem_correction(70, DMSO=3, DMSOfactor=0.6), 70 - 0.6 * 3
        )

    def test_formamide_method_one(self):
        """McConaughy et al. (1969): Tm drops by a factor per percent."""
        self.assertAlmostEqual(mt.chem_correction(70, fmd=5), 70 - 0.65 * 5)
        self.assertAlmostEqual(
            mt.chem_correction(70, fmd=5, fmdfactor=0.72), 70 - 0.72 * 5
        )

    def test_formamide_method_two(self):
        """Blake & Delcourt (1996): Tm += (0.453 f(GC) - 2.88) [formamide]."""
        for gc in (0, 25, 50, 100):
            with self.subTest(GC=gc):
                self.assertAlmostEqual(
                    mt.chem_correction(70, fmd=1.25, fmdmethod=2, GC=gc),
                    70 + (0.453 * (gc / 100.0) - 2.88) * 1.25,
                    places=10,
                )

    def test_formamide_method_two_needs_gc(self):
        """Method 2 needs a GC content."""
        with self.assertRaises(ValueError) as context:
            mt.chem_correction(70, fmd=1.25, fmdmethod=2)
        self.assertIn("GC", str(context.exception))
        with self.assertRaises(ValueError):
            mt.chem_correction(70, fmd=1.25, fmdmethod=2, GC=-1)

    def test_unknown_formamide_method(self):
        """Only formamide methods 1 and 2 exist."""
        with self.assertRaises(ValueError) as context:
            mt.chem_correction(70, fmd=1.25, fmdmethod=3, GC=50)
        self.assertIn("fmdmethod", str(context.exception))


class TmWallaceTests(unittest.TestCase):
    """Tests for Bio.SeqUtils.MeltingTemp.Tm_Wallace."""

    def test_wallace_rule(self):
        """Tm = 2 degC per A or T plus 4 degC per G or C."""
        for seq in ("ACGTTGCAATGCCGTA", "AAAAAAAA", "GCGCGCGC", "ATGC"):
            with self.subTest(seq=seq):
                expected = 2 * (seq.count("A") + seq.count("T")) + 4 * (
                    seq.count("G") + seq.count("C")
                )
                self.assertAlmostEqual(mt.Tm_Wallace(seq), expected)

    def test_strong_and_weak_are_not_ambiguous(self):
        """W counts as A/T and S counts as G/C, even when strict."""
        self.assertAlmostEqual(mt.Tm_Wallace("WWWW"), 8)
        self.assertAlmostEqual(mt.Tm_Wallace("SSSS"), 16)

    def test_ambiguous_bases_rejected_when_strict(self):
        """Truly ambiguous bases are refused by default."""
        for base in "BDHKMNRVY":
            with self.subTest(base=base):
                with self.assertRaises(ValueError):
                    mt.Tm_Wallace("ACGT" + base)

    def test_ambiguous_bases_averaged_when_not_strict(self):
        """Ambiguous bases contribute the average of what they stand for."""
        # N, K, M, R and Y are half A/T and half G/C: (2 + 4) / 2 = 3
        for base in "KMNRY":
            with self.subTest(base=base):
                self.assertAlmostEqual(
                    mt.Tm_Wallace("ACGT" + base, strict=False), 12 + 3
                )
        # B and V are two thirds G/C: (2 + 4 + 4) / 3
        for base in "BV":
            with self.subTest(base=base):
                self.assertAlmostEqual(
                    mt.Tm_Wallace("ACGT" + base, strict=False), 12 + 10 / 3.0
                )
        # D and H are two thirds A/T: (2 + 2 + 4) / 3
        for base in "DH":
            with self.subTest(base=base):
                self.assertAlmostEqual(
                    mt.Tm_Wallace("ACGT" + base, strict=False), 12 + 8 / 3.0
                )

    def test_whitespace_and_junk_ignored(self):
        """Non-base characters and whitespace do not count."""
        self.assertAlmostEqual(mt.Tm_Wallace("ACGT TGCA ATGC CGTA"), 48.0)
        self.assertAlmostEqual(mt.Tm_Wallace("1ACGT2TGCA3ATGC4CGTA"), 48.0)


class TmGCTests(unittest.TestCase):
    """Tests for Bio.SeqUtils.MeltingTemp.Tm_GC."""

    # (A, B, C, D) and the salt correction each value set implies, from the
    # references listed in the Tm_GC docstring.
    valuesets = {
        1: ((69.3, 0.41, 650, 1), 0),
        2: ((81.5, 0.41, 675, 1), 0),
        3: ((81.5, 0.41, 675, 1), 1),
        4: ((81.5, 0.41, 500, 1), 2),
        5: ((78.0, 0.7, 500, 1), 2),
        6: ((67.0, 0.8, 500, 1), 2),
        7: ((81.5, 0.41, 600, 1), 1),
        8: ((77.1, 0.41, 528, 1), 4),
    }
    seq = "CGTTCCAAAGATGTGGGCATGAGCTTAC"

    def reference(self, seq, valueset, Na=50, mismatch=True):
        """Return Tm = A + B(%GC) - C/N (+ salt) - D(%mismatch)."""
        (a, b, c, d), saltcorr = self.valuesets[valueset]
        percent_gc = 100 * (seq.count("G") + seq.count("C")) / len(seq)
        percent_gc += 50.0 * seq.count("X") / len(seq)  # gc_fraction counts X as 0.5
        if mismatch:
            percent_gc -= seq.count("X") * 50.0 / len(seq)
        melting_temp = a + b * percent_gc - c / len(seq)
        if saltcorr:
            melting_temp += reference_salt_correction(Na=Na, method=saltcorr, seq=seq)
        if mismatch:
            melting_temp -= d * (seq.count("X") * 100.0 / len(seq))
        return melting_temp

    def test_all_value_sets(self):
        """Each value set reproduces its published formula."""
        for valueset in self.valuesets:
            with self.subTest(valueset=valueset):
                self.assertAlmostEqual(
                    mt.Tm_GC(self.seq, valueset=valueset),
                    self.reference(self.seq, valueset),
                    places=10,
                )

    def test_value_set_overrides_salt_correction_argument(self):
        """A value set fixes its own salt correction method."""
        # value set 2 implies no salt correction at all
        self.assertAlmostEqual(
            mt.Tm_GC(self.seq, valueset=2, saltcorr=4),
            mt.Tm_GC(self.seq, valueset=2, saltcorr=0),
            places=10,
        )
        # ...and the result is salt independent
        self.assertAlmostEqual(
            mt.Tm_GC(self.seq, valueset=2, Na=50),
            mt.Tm_GC(self.seq, valueset=2, Na=500),
            places=10,
        )

    def test_userset_overrides_value_set(self):
        """A user supplied (A, B, C, D) is used instead of a value set."""
        userset = (70.0, 0.5, 600, 2)
        percent_gc = 100 * (self.seq.count("G") + self.seq.count("C")) / len(self.seq)
        expected = 70.0 + 0.5 * percent_gc - 600 / len(self.seq)
        self.assertAlmostEqual(
            mt.Tm_GC(self.seq, userset=userset, saltcorr=0), expected, places=10
        )

    def test_mismatch_penalty(self):
        """Every X costs D degrees per percent of the sequence."""
        seq = "CGTTCCAAAGATGTXGGCATGAGCTTAC"
        self.assertAlmostEqual(
            mt.Tm_GC(seq, valueset=7, strict=False),
            self.reference(seq, 7, mismatch=True),
            places=10,
        )
        # with mismatch=False the X still counts as half a G/C
        (a, b, c, _), saltcorr = self.valuesets[7]
        percent_gc = 100 * (seq.count("G") + seq.count("C") + 0.5) / len(seq)
        expected = (
            a
            + b * percent_gc
            - c / len(seq)
            + reference_salt_correction(Na=50, method=saltcorr, seq=seq)
        )
        self.assertAlmostEqual(
            mt.Tm_GC(seq, valueset=7, strict=False, mismatch=False),
            expected,
            places=10,
        )

    def test_salt_correction_five_rejected(self):
        """Method 5 corrects entropy, so it cannot correct a Tm directly."""
        with self.assertRaises(ValueError) as context:
            mt.Tm_GC(self.seq, saltcorr=5)
        self.assertIn("5", str(context.exception))

    def test_ambiguous_bases_rejected_when_strict(self):
        """Ambiguous bases are refused by default."""
        for base in "BDHKMNRVY":
            with self.subTest(base=base):
                with self.assertRaises(ValueError):
                    mt.Tm_GC("CGTTCCAAAGATGTGGGCATGAGCTTA" + base)

    def test_unknown_value_set_rejected(self):
        """Value sets above 8 do not exist."""
        with self.assertRaises(ValueError) as context:
            mt.Tm_GC(self.seq, valueset=9)
        self.assertIn("valueset", str(context.exception))

    def test_tm_rises_with_gc_content(self):
        """More G+C means a higher predicted Tm."""
        temperatures = [
            mt.Tm_GC("A" * (20 - n) + "G" * n, strict=False) for n in range(0, 21, 5)
        ]
        self.assertEqual(temperatures, sorted(temperatures))


class TmNNTests(unittest.TestCase):
    """Tests for Bio.SeqUtils.MeltingTemp.Tm_NN."""

    sequences = [
        "CGTTCCAAAGATGTGGGCATGAGCTTAC",
        "ACGTTGCAATGCCGTA",
        "TTTTTTTTTTAAAA",  # no G or C at all, starts with T
        "AAAAAAAAAAAA",  # no G or C at all, ends with A
        "TGCAGCTAGCTA",  # starts with T and ends with A
        "GGGGCCCCGGGG",
        "ATGC",
    ]

    def test_perfect_match_against_reference(self):
        """A perfectly matched duplex agrees with the NN model summation."""
        for table_name in ("DNA_NN1", "DNA_NN2", "DNA_NN3", "DNA_NN4"):
            table = getattr(mt, table_name)
            for seq in self.sequences:
                with self.subTest(table=table_name, seq=seq):
                    self.assertAlmostEqual(
                        mt.Tm_NN(seq, nn_table=table, saltcorr=0),
                        reference_tm_nn(seq, nn_table=table),
                        places=8,
                    )

    def test_rna_and_hybrid_tables(self):
        """The RNA and RNA/DNA tables are summed the same way."""
        for table_name in ("RNA_NN1", "RNA_NN2", "RNA_NN3", "R_DNA_NN1"):
            table = getattr(mt, table_name)
            for seq in ("ACGTTGCAATGCCGTA", "CGTTCCAAAGATGTGGGCATGAGCTTAC"):
                with self.subTest(table=table_name, seq=seq):
                    self.assertAlmostEqual(
                        mt.Tm_NN(seq, nn_table=table, saltcorr=0),
                        reference_tm_nn(seq, nn_table=table),
                        places=8,
                    )

    def test_concentrations(self):
        """The strand concentration enters as R ln(dnac1 - dnac2/2)."""
        seq = "ACGTTGCAATGCCGTA"
        for dnac1, dnac2 in ((25, 25), (100, 50), (500, 0), (50, 10)):
            with self.subTest(dnac1=dnac1, dnac2=dnac2):
                self.assertAlmostEqual(
                    mt.Tm_NN(seq, dnac1=dnac1, dnac2=dnac2, saltcorr=0),
                    reference_tm_nn(seq, dnac1=dnac1, dnac2=dnac2),
                    places=8,
                )
        self.assertLess(
            mt.Tm_NN(seq, dnac1=25, dnac2=25, saltcorr=0),
            mt.Tm_NN(seq, dnac1=250, dnac2=250, saltcorr=0),
        )

    def test_self_complementary(self):
        """A self complementary duplex uses dnac1 and the symmetry term."""
        seq = "ACGCGCGT"  # its own reverse complement
        self.assertEqual(seq, reverse_complement(seq))
        self.assertAlmostEqual(
            mt.Tm_NN(seq, selfcomp=True, saltcorr=0),
            reference_tm_nn(seq, selfcomp=True),
            places=8,
        )
        # the symmetry correction lowers the entropy, so it raises Tm here
        self.assertNotAlmostEqual(
            mt.Tm_NN(seq, selfcomp=True, saltcorr=0),
            mt.Tm_NN(seq, selfcomp=False, saltcorr=0),
        )

    def test_strand_symmetry(self):
        """A duplex melts at the same temperature from either strand."""
        for seq in self.sequences:
            with self.subTest(seq=seq):
                self.assertAlmostEqual(
                    mt.Tm_NN(seq, saltcorr=0),
                    mt.Tm_NN(reverse_complement(seq), saltcorr=0),
                    places=8,
                )

    def test_salt_corrections_are_applied_in_the_right_place(self):
        """Methods 1-4 shift Tm, 5 shifts dS and 6-7 shift 1/Tm."""
        seq = "CGTTCCAAAGATGTGGGCATGAGCTTAC"
        uncorrected = reference_tm_nn(seq)
        for method in (1, 2, 3, 4):
            with self.subTest(method=method):
                correction = reference_salt_correction(Na=50, method=method, seq=seq)
                self.assertAlmostEqual(
                    mt.Tm_NN(seq, saltcorr=method),
                    uncorrected + correction,
                    places=8,
                )
        correction = reference_salt_correction(Na=50, method=5, seq=seq)
        self.assertAlmostEqual(
            mt.Tm_NN(seq, saltcorr=5),
            reference_tm_nn(seq, extra=(0.0, correction)),
            places=8,
        )
        for method in (6, 7):
            with self.subTest(method=method):
                correction = reference_salt_correction(Na=50, method=method, seq=seq)
                self.assertAlmostEqual(
                    mt.Tm_NN(seq, saltcorr=method),
                    1 / (1 / (uncorrected + 273.15) + correction) - 273.15,
                    places=8,
                )

    def test_internal_mismatch(self):
        """An internal mismatch is looked up in the internal mismatch table."""
        seq = "ATGGCCATTGTAA"
        c_seq = "TACCGATAACATT"  # complement, with position 5 mismatched
        self.assertEqual(len(seq), len(c_seq))
        self.assertAlmostEqual(
            mt.Tm_NN(seq, c_seq=c_seq, saltcorr=0),
            reference_tm_nn(seq, c_seq=c_seq, imm_table=mt.DNA_IMM1),
            places=8,
        )
        # a mismatch must destabilise the duplex
        self.assertLess(
            mt.Tm_NN(seq, c_seq=c_seq, saltcorr=0), mt.Tm_NN(seq, saltcorr=0)
        )

    def test_left_terminal_mismatch(self):
        """A mismatch at the 5' end uses the terminal mismatch table."""
        seq = "ATGGCCATTGTAA"
        c_seq = "AACCGGTAACATT"  # complement with the first base mismatched
        expected = reference_tm_nn(
            seq[1:],
            c_seq[1:],
            init_seq=seq,
            extra=mt.DNA_TMM1["AA/TA"],
        )
        self.assertAlmostEqual(
            mt.Tm_NN(seq, c_seq=c_seq, saltcorr=0), expected, places=8
        )

    def test_right_terminal_mismatch(self):
        """A mismatch at the 3' end uses the terminal mismatch table."""
        seq = "ATGGCCATTGTAA"
        c_seq = "TACCGGTAACATA"  # complement with the last base mismatched
        expected = reference_tm_nn(
            seq[:-1],
            c_seq[:-1],
            init_seq=seq,
            extra=mt.DNA_TMM1["AA/TA"],
        )
        self.assertAlmostEqual(
            mt.Tm_NN(seq, c_seq=c_seq, saltcorr=0), expected, places=8
        )

    def test_dangling_end_on_the_template(self):
        """A template overhanging both ends adds two dangling end terms."""
        seq = "ACGTTGCA"
        c_seq = "G" + "TGCAACGT" + "C"
        extra = (
            mt.DNA_DE1[".A/GT"][0] + mt.DNA_DE1["CT/.A"][0],
            mt.DNA_DE1[".A/GT"][1] + mt.DNA_DE1["CT/.A"][1],
        )
        expected = reference_tm_nn(seq, "TGCAACGT", init_seq=seq, extra=extra)
        self.assertAlmostEqual(
            mt.Tm_NN(seq, c_seq=c_seq, shift=1, saltcorr=0), expected, places=8
        )

    def test_negative_shift(self):
        """A primer overhanging both ends adds two dangling end terms."""
        seq = "GACGTTGCAC"
        c_seq = "TGCAACGT"
        extra = (
            mt.DNA_DE1["GA/.T"][0] + mt.DNA_DE1[".T/CA"][0],
            mt.DNA_DE1["GA/.T"][1] + mt.DNA_DE1[".T/CA"][1],
        )
        expected = reference_tm_nn("ACGTTGCA", "TGCAACGT", init_seq=seq, extra=extra)
        self.assertAlmostEqual(
            mt.Tm_NN(seq, c_seq=c_seq, shift=-1, saltcorr=0), expected, places=8
        )

    def test_short_template_gives_a_single_dangling_end(self):
        """A template shorter at the 3' end leaves one dangling base."""
        seq = "ACGTTGCAAG"
        c_seq = "TGCAACGT"  # complement of seq[:8] only
        expected = reference_tm_nn(
            "ACGTTGCA", "TGCAACGT", init_seq=seq, extra=mt.DNA_DE1[".T/AA"]
        )
        self.assertAlmostEqual(
            mt.Tm_NN(seq, c_seq=c_seq, saltcorr=0), expected, places=8
        )

    def test_over_dangling_ends_are_trimmed(self):
        """Only one dangling base per end is kept."""
        seq = "ACGTTGCA"
        c_seq = "GT" + "TGCAACGT"  # two extra bases, of which one is dropped
        expected = reference_tm_nn(
            seq, "TGCAACGT", init_seq=seq, extra=mt.DNA_DE1[".A/TT"]
        )
        self.assertAlmostEqual(
            mt.Tm_NN(seq, c_seq=c_seq, shift=2, saltcorr=0), expected, places=8
        )

    def test_missing_dangling_end_data_is_reported(self):
        """A dangling end with no thermodynamic data is reported."""
        seq = "ACGTTGCAAG"
        c_seq = "TGCAACGA"  # last base does not pair, so ".A/AA" is needed
        with self.assertRaises(ValueError) as context:
            mt.Tm_NN(seq, c_seq=c_seq, saltcorr=0)
        self.assertIn(".A/AA", str(context.exception))

    def test_missing_neighbour_data_is_an_error_when_strict(self):
        """Unknown neighbours raise by default."""
        with self.assertRaises(ValueError) as context:
            mt.Tm_NN("ACGTACGT", c_seq="TGCANGCA", check=False)
        self.assertIn("no thermodynamic data", str(context.exception))

    def test_missing_neighbour_data_only_warns_when_not_strict(self):
        """Unknown neighbours only warn when strict is False."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            value = mt.Tm_NN("ACGTACGT", c_seq="TGCANGCA", check=False, strict=False)
        self.assertIsInstance(value, float)
        self.assertTrue(
            any(issubclass(entry.category, BiopythonWarning) for entry in caught)
        )

    def test_tm_rises_with_gc_content(self):
        """More G+C means a higher predicted Tm."""
        temperatures = [
            mt.Tm_NN("AT" * (8 - n) + "GC" * n, saltcorr=0) for n in range(9)
        ]
        self.assertEqual(temperatures, sorted(temperatures))

    def test_accepts_seq_objects_and_whitespace(self):
        """Seq objects, RNA, whitespace and lower case are all accepted."""
        expected = mt.Tm_NN("ACGTTGCAATGCCGTA", saltcorr=0)
        for value in (
            Seq("ACGTTGCAATGCCGTA"),
            "acgttgcaatgccgta",
            "ACGT TGCA ATGC CGTA",
            "ACGUUGCAAUGCCGUA",
        ):
            with self.subTest(value=str(value)):
                self.assertAlmostEqual(mt.Tm_NN(value, saltcorr=0), expected, places=8)


class MakeTableTests(unittest.TestCase):
    """Tests for Bio.SeqUtils.MeltingTemp.make_table."""

    def test_default_table_is_all_zeroes(self):
        """Without an old table every entry starts at (0, 0)."""
        table = mt.make_table()
        self.assertEqual(set(table.values()), {(0, 0)})
        self.assertIn("init", table)
        self.assertIn("GG/CC", table)

    def test_values_replace_entries(self):
        """The values argument replaces individual entries."""
        table = mt.make_table(values={"init_A/T": (2.3, 4.1)})
        self.assertEqual(table["init_A/T"], (2.3, 4.1))
        self.assertEqual(table["init_G/C"], (0, 0))

    def test_old_table_is_not_modified(self):
        """Updating a table leaves the original untouched."""
        before = dict(mt.DNA_NN2)
        table = mt.make_table(oldtable=mt.DNA_NN2, values={"init_A/T": (2.3, 4.1)})
        self.assertEqual(mt.DNA_NN2, before)
        self.assertEqual(table["init_A/T"], (2.3, 4.1))
        self.assertEqual(table["AA/TT"], mt.DNA_NN2["AA/TT"])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
