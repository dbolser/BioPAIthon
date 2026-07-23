# Copyright 2003 by Iddo Friedberg.  All rights reserved.
# Copyright 2007-2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for SeqUtils module."""

import os
import re
import unittest

from Bio import SeqIO
from Bio.Data import IUPACData
from Bio.Seq import complement
from Bio.Seq import complement_rna
from Bio.Seq import MutableSeq
from Bio.Seq import reverse_complement
from Bio.Seq import Seq
from Bio.Seq import translate
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import CodonAdaptationIndex
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import GC123
from Bio.SeqUtils import GC_skew
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils import nt_search
from Bio.SeqUtils import seq1
from Bio.SeqUtils import seq3
from Bio.SeqUtils import six_frame_translations
from Bio.SeqUtils.CheckSum import crc32
from Bio.SeqUtils.CheckSum import crc64
from Bio.SeqUtils.CheckSum import gcg
from Bio.SeqUtils.CheckSum import seguid
from Bio.SeqUtils.lcc import lcc_mult
from Bio.SeqUtils.lcc import lcc_simp


class SeqUtilsTests(unittest.TestCase):
    # Example of crc64 collision from Sebastian Bassi using the
    # immunoglobulin lambda light chain variable region from Homo sapiens
    # Both sequences share the same CRC64 checksum: 44CAAD88706CC153
    str_light_chain_one = (
        "QSALTQPASVSGSPGQSITISCTGTSSDVGSYNLVSWYQQHPGKAPKLMIYEGSKRPSGV"
        "SNRFSGSKSGNTASLTISGLQAEDEADYYCSSYAGSSTLVFGGGTKLTVL"
    )
    str_light_chain_two = (
        "QSALTQPASVSGSPGQSITISCTGTSSDVGSYNLVSWYQQHPGKAPKLMIYEGSKRPSGV"
        "SNRFSGSKSGNTASLTISGLQAEDEADYYCCSYAGSSTWVFGGGTKLTVL"
    )

    def test_codon_adaptation_index_initialization(self):
        """Test Codon Adaptation Index (CAI) initialization from sequences."""
        # We need CDS sequences to count the codon usage...
        dna_filename = "GenBank/NC_005816.gb"
        record = SeqIO.read(dna_filename, "genbank")
        records = []
        for feature in record.features:
            if feature.type == "CDS" and len(feature.location.parts) == 1:
                start = feature.location.start
                end = feature.location.end
                table = int(feature.qualifiers["transl_table"][0])
                if feature.location.strand == -1:
                    seq = record.seq[start:end].reverse_complement()
                else:
                    seq = record.seq[start:end]
                # Double check we have the CDS sequence expected
                # TODO - Use any cds_start option if/when added to deal with the met
                a = "M" + seq[3:].translate(table)
                b = feature.qualifiers["translation"][0] + "*"
                self.assertEqual(a, b)
                records.append(
                    SeqRecord(
                        seq,
                        id=feature.qualifiers["protein_id"][0],
                        description=feature.qualifiers["product"][0],
                    )
                )

        cai = CodonAdaptationIndex(records)
        # Now check codon usage index (CAI) using this species
        self.assertEqual(
            record.annotations["source"], "Yersinia pestis biovar Microtus str. 91001"
        )
        value = cai.calculate("ATGCGTATCGATCGCGATACGATTAGGCGGATG")
        self.assertAlmostEqual(value, 0.70246, places=5)
        optimized_sequence = cai.optimize(
            "ATGCGTATCGATCGCGATACGATTAGGCGGATG", strict=False
        )
        optimized_value = cai.calculate(optimized_sequence)
        self.assertEqual(optimized_value, 1.0)
        aa_initial = Seq("ATGCGTATCGATCGCGATACGATTAGGCGGATG").translate()
        aa_optimized = optimized_sequence.translate()
        self.assertEqual(aa_initial, aa_optimized)
        with self.assertRaises(KeyError):
            cai.optimize("CAU", "protein", strict=False)
        self.assertEqual(
            str(cai),
            """\
AAA	1.000
AAC	0.385
AAG	0.344
AAT	1.000
ACA	1.000
ACC	0.553
ACG	0.319
ACT	0.447
AGA	0.595
AGC	0.967
AGG	0.297
AGT	1.000
ATA	0.581
ATC	0.930
ATG	1.000
ATT	1.000
CAA	0.381
CAC	0.581
CAG	1.000
CAT	1.000
CCA	0.500
CCC	0.500
CCG	1.000
CCT	0.767
CGA	0.568
CGC	0.919
CGG	0.514
CGT	1.000
CTA	0.106
CTC	0.379
CTG	1.000
CTT	0.424
GAA	1.000
GAC	0.633
GAG	0.506
GAT	1.000
GCA	1.000
GCC	0.617
GCG	0.532
GCT	0.809
GGA	1.000
GGC	0.525
GGG	0.575
GGT	0.950
GTA	0.500
GTC	0.618
GTG	0.971
GTT	1.000
TAA	1.000
TAC	0.434
TAG	0.062
TAT	1.000
TCA	1.000
TCC	0.533
TCG	0.233
TCT	0.967
TGA	0.250
TGC	1.000
TGG	1.000
TGT	0.750
TTA	0.455
TTC	1.000
TTG	0.212
TTT	0.886
""",
        )

    def test_codon_adaptation_index_calculation(self):
        """Test Codon Adaptation Index (CAI) calculation for an mRNA."""
        cai = CodonAdaptationIndex([])
        # Use the Codon Adaption Index for E. coli, precalculated by
        # Sharp and Li (Nucleic Acids Res. 1987 Feb 11;15(3):1281-95), Table 1.
        cai["TTT"] = 0.296  # Phe
        cai["TTC"] = 1.000  # Phe
        cai["TTA"] = 0.020  # Leu
        cai["TTG"] = 0.020  # Leu
        cai["CTT"] = 0.042  # Leu
        cai["CTC"] = 0.037  # Leu
        cai["CTA"] = 0.007  # Leu
        cai["CTG"] = 1.000  # Leu
        cai["ATT"] = 0.185  # Ile
        cai["ATC"] = 1.000  # Ile
        cai["ATA"] = 0.003  # Ile
        cai["ATG"] = 1.000  # Met
        cai["GTT"] = 1.000  # Val
        cai["GTC"] = 0.066  # Val
        cai["GTA"] = 0.495  # Val
        cai["GTG"] = 0.221  # Val
        cai["TAT"] = 0.239  # Tyr
        cai["TAC"] = 1.000  # Tyr
        cai["CAT"] = 0.291  # His
        cai["CAC"] = 1.000  # His
        cai["CAA"] = 0.124  # Gln
        cai["CAG"] = 1.000  # Gln
        cai["AAT"] = 0.051  # Asn
        cai["AAC"] = 1.000  # Asn
        cai["AAA"] = 1.000  # Lys
        cai["AAG"] = 0.253  # Lys
        cai["GAT"] = 0.434  # Asp
        cai["GAC"] = 1.000  # Asp
        cai["GAA"] = 1.000  # Glu
        cai["GAG"] = 0.259  # Glu
        cai["TCT"] = 1.000  # Ser
        cai["TCC"] = 0.744  # Ser
        cai["TCA"] = 0.077  # Ser
        cai["TCG"] = 0.017  # Ser
        cai["CCT"] = 0.070  # Pro
        cai["CCC"] = 0.012  # Pro
        cai["CCA"] = 0.135  # Pro
        cai["CCG"] = 1.000  # Pro
        cai["ACT"] = 0.965  # Thr
        cai["ACC"] = 1.000  # Thr
        cai["ACA"] = 0.076  # Thr
        cai["ACG"] = 0.099  # Thr
        cai["GCT"] = 1.000  # Ala
        cai["GCC"] = 0.122  # Ala
        cai["GCA"] = 0.586  # Ala
        cai["GCG"] = 0.424  # Ala
        cai["TGT"] = 0.500  # Cys
        cai["TGC"] = 1.000  # Cys
        cai["TGG"] = 1.000  # Trp
        cai["CGT"] = 1.000  # Arg
        cai["CGC"] = 0.356  # Arg
        cai["CGA"] = 0.004  # Arg
        cai["CGG"] = 0.004  # Arg
        cai["AGT"] = 0.085  # Ser
        cai["AGC"] = 0.410  # Ser
        cai["AGA"] = 0.004  # Arg
        cai["AGG"] = 0.002  # Arg
        cai["GGT"] = 1.000  # Gly
        cai["GGC"] = 0.724  # Gly
        cai["GGA"] = 0.010  # Gly
        cai["GGG"] = 0.019  # Gly
        # Now calculate the CAI for the genes listed in Table 2 of
        # Sharp and Li (Nucleic Acids Res. 1987 Feb 11;15(3):1281-95).
        rpsU = Seq(
            "CCGGTAATTAAAGTACGTGAAAACGAGCCGTTCGACGTAGCTCTGCGTCGCTTCAAGCGTTCCTGCGAAAAAGCAGGTGTTCTGGCGGAAGTTCGTCGTCGTGAGTTCTATGAAAAACCGACTACCGAACGTAAGCGCGCTAAAGCTTCTGCAGTGAAACGTCACGCGAAGAAACTGGCTCGCGAAAACGCACGCCGCACTCGTCTGTAC"
        )
        self.assertAlmostEqual(cai.calculate(rpsU), 0.726, places=3)
        rpoD = Seq(
            "ATGGAGCAAAACCCGCAGTCACAGCTGAAACTTCTTGTCACCCGTGGTAAGGAGCAAGGCTATCTGACCTATGCCGAGGTCAATGACCATCTGCCGGAAGATATCGTCGATTCAGATCAGATCGAAGACATCATCCAAATGATCAACGACATGGGCATTCAGGTGATGGAAGAAGCACCGGATGCCGATGATCTGATGCTGGCTGAAAACACCGCGGACGAAGATGCTGCCGAAGCCGCCGCGCAGGTGCTTTCCAGCGTGGAATCTGAAATCGGGCGCACGACTGACCCGGTACGCATGTACATGCGTGAAATGGGCACCGTTGAACTGTTGACCCGCGAAGGCGAAATTGACATCGCTAAGCGTATTGAAGACGGGATCAACCAGGTTCAATGCTCCGTTGCTGAATATCCGGAAGCGATCACCTATCTGCTGGAACAGTACGATCGTGTTGAAGCAGAAGAAGCGCGTCTGTCCGATCTGATCACCGGCTTTGTTGACCCGAACGCAGAAGAAGATCTGGCACCTACCGCCACTCACGTCGGTTCTGAGCTTTCCCAGGAAGATCTGGACGATGACGAAGATGAAGACGAAGAAGATGGCGATGACGACAGCGCCGATGATGACAACAGCATCGACCCGGAACTGGCTCGCGAAAAATTTGCGGAACTACGCGCTCAGTACGTTGTAACGCGTGACACCATCAAAGCGAAAGGTCGCAGTCACGCTACCGCTCAGGAAGAGATCCTGAAACTGTCTGAAGTATTCAAACAGTTCCGCCTGGTGCCGAAGCAGTTTGACTACCTGGTCAACAGCATGCGCGTCATGATGGACCGCGTTCGTACGCAAGAACGTCTGATCATGAAGCTCTGCGTTGAGCAGTGCAAAATGCCGAAGAAAAACTTCATTACCCTGTTTACCGGCAACGAAACCAGCGATACCTGGTTCAACGCGGCAATTGCGATGAACAAGCCGTGGTCGGAAAAACTGCACGATGTCTCTGAAGAAGTGCATCGCGCCCTGCAAAAACTGCAGCAGATTGAAGAAGAAACCGGCCTGACCATCGAGCAGGTTAAAGATATCAACCGTCGTATGTCCATCGGTGAAGCGAAAGCCCGCCGTGCGAAGAAAGAGATGGTTGAAGCGAACTTACGTCTGGTTATTTCTATCGCTAAGAAATACACCAACCGTGGCTTGCAGTTCCTTGACCTGATTCAGGAAGGCAACATCGGTCTGATGAAAGCGGTTGATAAATTCGAATACCGCCGTGGTTACAAGTTCTCCACCTACGCAACCTGGTGGATCCGTCAGGCGATCACCCGCTCTATCGCGGATCAGGCGCGCACCATCCGTATTCCGGTGCATATGATTGAGACCATCAACAAGCTCAACCGTATTTCTCGCCAGATGCTGCAAGAGATGGGCCGTGAACCGACGCCGGAAGAACTGGCTGAACGTATGCTGATGCCGGAAGACAAGATCCGCAAAGTGCTGAAGATCGCCAAAGAGCCAATCTCCATGGAAACGCCGATCGGTGATGATGAAGATTCGCATCTGGGGGATTTCATCGAGGATACCACCCTCGAGCTGCCGCTGGATTCTGCGACCACCGAAAGCCTGCGTGCGGCAACGCACGACGTGCTGGCTGGCCTGACCGCGCGTGAAGCAAAAGTTCTGCGTATGCGTTTCGGTATCGATATGAACACCGACTACACGCTGGAAGAAGTGGGTAAACAGTTCGACGTTACCCGCGAACGTATCCGTCAGATCGAAGCGAAGGCGCTGCGCAAACTGCGTCACCCGAGCCGTTCTGAAGTGCTGCGTAGCTTCCTGGACGAT"
        )
        self.assertAlmostEqual(cai.calculate(rpoD), 0.582, places=2)
        dnaG = "ATGGCTGGACGAATCCCACGCGTATTCATTAATGATCTGCTGGCACGCACTGACATCGTCGATCTGATCGATGCCCGTGTGAAGCTGAAAAAGCAGGGCAAGAATTTCCACGCGTGTTGTCCATTCCACAACGAGAAAACCCCGTCCTTCACCGTTAACGGTGAGAAACAGTTTTACCACTGCTTTGGATGTGGCGCGCACGGCAACGCGATCGACTTCCTGATGAACTACGACAAGCTCGAGTTCGTCGAAACGGTCGAAGAGCTGGCAGCAATGCACAATCTTGAAGTGCCATTTGAAGCAGGCAGCGGCCCCAGCCAGATCGAGCGCCATCAGAGGCAAACGCTTTATCAGTTGATGGACGGTCTGAATACGTTTTACCAACAATCTTTACAACAACCTGTTGCCACGTCTGCGCGCCAGTATCTGGAAAAACGCGGATTAAGCCACGAGGTTATCGCTCGCTTTGCGATTGGTTTTGCGCCCCCCGGCTGGGACAACGTCCTGAAGCGGTTTGGCGGCAATCCAGAAAATCGCCAGTCATTGATTGATGCGGGGATGTTGGTCACTAACGATCAGGGACGCAGTTACGATCGTTTCCGCGAGCGGGTGATGTTCCCCATTCGCGATAAACGCGGTCGGGTGATTGGTTTTGGCGGGCGCGTGCTGGGCAACGATACCCCCAAATACCTGAACTCGCCGGAAACAGACATTTTCCATAAAGGCCGCCAGCTTTACGGTCTTTATGAAGCGCAGCAGGATAACGCTGAACCCAATCGTCTGCTTGTGGTCGAAGGCTATATGGACGTGGTGGCGCTGGCGCAATACGGCATTAATTACGCCGTTGCGTCGTTAGGTACGTCAACCACCGCCGATCACATACAACTGTTGTTCCGCGCGACCAACAATGTCATTTGCTGTTATGACGGCGACCGTGCAGGCCGCGATGCCGCCTGGCGAGCGCTGGAAACGGCGCTGCCTTACATGACAGACGGCCGTCAGCTACGCTTTATGTTTTTGCCTGATGGCGAAGACCCTGACACGCTAGTACGAAAAGAAGGTAAAGAAGCGTTTGAAGCGCGGATGGAGCAGGCGATGCCACTCTCCGCATTTCTGTTTAACAGTCTGATGCCGCAAGTTGATCTGAGTACCCCTGACGGGCGCGCACGTTTGAGTACGCTGGCACTACCATTGATATCGCAAGTGCCGGGCGAAACGCTGCGAATATATCTTCGTCAGGAATTAGGCAACAAATTAGGCATACTTGATGACAGCCAGCTTGAACGATTAATGCCAAAAGCGGCAGAGAGCGGCGTTTCTCGCCCTGTTCCGCAGCTAAAACGCACGACCATGCGTATACTTATAGGGTTGCTGGTGCAAAATCCAGAATTAGCGACGTTGGTCCCGCCGCTTGAGAATCTGGATGAAAATAAGCTCCCTGGACTTGGCTTATTCAGAGAACTGGTCAACACTTGTCTCTCCCAGCCAGGTCTGACCACCGGGCAACTTTTAGAGCACTATCGTGGTACAAATAATGCTGCCACCCTTGAAAAACTGTCGATGTGGGACGATATAGCAGATAAGAATATTGCTGAGCAAACCTTCACCGACTCACTCAACCATATGTTTGATTCGCTGCTTGAACTGCGCCAGGAAGAGTTAATCGCTCGTGAGCGCACGCATGGTTTAAGCAACGAAGAACGCCTGGAGCTCTGGACATTAAACCAGGAGCTGGCGAAAAAG"
        self.assertAlmostEqual(cai.calculate(dnaG), 0.271, places=3)
        lacI = "GTGAAACCAGTAACGTTATACGATGTCGCAGAGTATGCCGGTGTCTCTTATCAGACCGTTTCCCGCGTGGTGAACCAGGCCAGCCACGTTTCTGCGAAAACGCGGGAAAAAGTGGAAGCGGCGATGGCGGAGCTGAATTACATTCCCAACCGCGTGGCACAACAACTGGCGGGCAAACAGTCGTTGCTGATTGGCGTTGCCACCTCCAGTCTGGCCCTGCACGCGCCGTCGCAAATTGTCGCGGCGATTAAATCTCGCGCCGATCAACTGGGTGCCAGCGTGGTGGTGTCGATGGTAGAACGAAGCGGCGTCGAAGCCTGTAAAGCGGCGGTGCACAATCTTCTCGCGCAACGCGTCAGTGGGCTGATCATTAACTATCCGCTGGATGACCAGGATGCCATTGCTGTGGAAGCTGCCTGCACTAATGTTCCGGCGTTATTTCTTGATGTCTCTGACCAGACACCCATCAACAGTATTATTTTCTCCCATGAAGACGGTACGCGACTGGGCGTGGAGCATCTGGTCGCATTGGGTCACCAGCAAATCGCGCTGTTAGCGGGCCCATTAAGTTCTGTCTCGGCGCGTCTGCGTCTGGCTGGCTGGCATAAATATCTCACTCGCAATCAAATTCAGCCGATAGCGGAACGGGAAGGCGACTGGAGTGCCATGTCCGGTTTTCAACAAACCATGCAAATGCTGAATGAGGGCATCGTTCCCACTGCGATGCTGGTTGCCAACGATCAGATGGCGCTGGGCGCAATGCGCGCCATTACCGAGTCCGGGCTGCGCGTTGGTGCGGATATCTCGGTAGTGGGATACGACGATACCGAAGACAGCTCATGTTATATCCCGCCGTTAACCACCATCAAACAGGATTTTCGCCTGCTGGGGCAAACCAGCGTGGACCGCTTGCTGCAACTCTCTCAGGGCCAGGCGGTGAAGGGCAATCAGCTGTTGCCCGTCTCACTGGTGAAAAGAAAAACCACCCTGGCGCCCAATACGCAAACCGCCTCTCCCCGCGCGTTGGCCGATTCATTAATGCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCGGGCAG"
        self.assertAlmostEqual(cai.calculate(lacI), 0.296, places=2)
        trpR = "ATGGCCCAACAATCACCCTATTCAGCAGCGATGGCAGAACAGCGTCACCAGGAGTGGTTACGTTTTGTCGACCTGCTTAAGAATGCCTACCAAAACGATCTCCATTTACCGTTGTTAAACCTGATGCTGACGCCAGATGAGCGCGAAGCGTTGGGGACTCGCGTGCGTATTGTCGAAGAGCTGTTGCGCGGCGAAATGAGCCAGCGTGAGTTAAAAAATGAACTCGGCGCAGGCATCGCGACGATTACGCGTGGATCTAACAGCCTGAAAGCCGCGCCCGTCGAGCTGCGCCAGTGGCTGGAAGAGGTGTTGCTGAAAAGCGAT"
        self.assertAlmostEqual(cai.calculate(trpR), 0.267, places=2)
        lpp = "ATGAAAGCTACTAAACTGGTACTGGGCGCGGTAATCCTGGGTTCTACTCTGCTGGCAGGTTGCTCCAGCAACGCTAAAATCGATCAGCTGTCTTCTGACGTTCAGACTCTGAACGCTAAAGTTGACCAGCTGAGCAACGACGTGAACGCAATGCGTTCCGACGTTCAGGCTGCTAAAGATGACGCAGCTCGTGCTAACCAGCGTCTGGACAACATGGCTACTAAATACCGCAAG"
        self.assertAlmostEqual(cai.calculate(lpp), 0.849, places=3)

    def test_crc_checksum_collision(self):
        # Explicit testing of crc64 collision:
        self.assertNotEqual(self.str_light_chain_one, self.str_light_chain_two)
        self.assertNotEqual(
            crc32(self.str_light_chain_one), crc32(self.str_light_chain_two)
        )
        self.assertEqual(
            crc64(self.str_light_chain_one), crc64(self.str_light_chain_two)
        )
        self.assertNotEqual(
            gcg(self.str_light_chain_one), gcg(self.str_light_chain_two)
        )
        self.assertNotEqual(
            seguid(self.str_light_chain_one), seguid(self.str_light_chain_two)
        )

    def seq_checksums(
        self,
        seq_str,
        exp_crc32,
        exp_crc64,
        exp_gcg,
        exp_seguid,
        exp_simple_LCC,
        exp_window_LCC,
    ):
        for s in [seq_str, Seq(seq_str), MutableSeq(seq_str)]:
            self.assertEqual(exp_crc32, crc32(s))
            self.assertEqual(exp_crc64, crc64(s))
            self.assertEqual(exp_gcg, gcg(s))
            self.assertEqual(exp_seguid, seguid(s))
            self.assertAlmostEqual(exp_simple_LCC, lcc_simp(s), places=4)
            values = lcc_mult(s, 20)
            self.assertEqual(len(exp_window_LCC), len(values), values)
            for value1, value2 in zip(exp_window_LCC, values):
                self.assertAlmostEqual(value1, value2, places=2)

    def test_checksum1(self):
        self.seq_checksums(
            self.str_light_chain_one,
            2994980265,
            "CRC-44CAAD88706CC153",
            9729,
            "BpBeDdcNUYNsdk46JoJdw7Pd3BI",
            0.5160,
            (
                0.4982,
                0.4794,
                0.4794,
                0.4794,
                0.3241,
                0.2160,
                0.1764,
                0.1764,
                0.1764,
                0.1764,
                0.2657,
                0.2948,
                0.1287,
            ),
        )

    def test_checksum2(self):
        self.seq_checksums(
            self.str_light_chain_two,
            802105214,
            "CRC-44CAAD88706CC153",
            9647,
            "X5XEaayob1nZLOc7eVT9qyczarY",
            0.5343,
            (
                0.4982,
                0.4794,
                0.4794,
                0.4794,
                0.3241,
                0.2160,
                0.1764,
                0.1764,
                0.1764,
                0.1764,
                0.2657,
                0.2948,
                0.1287,
            ),
        )

    def test_checksum3(self):
        self.seq_checksums(
            "ATGCGTATCGATCGCGATACGATTAGGCGGAT",
            817679856,
            "CRC-6234FF451DC6DFC6",
            7959,
            "8WCUbVjBgiRmM10gfR7XJNjbwnE",
            0.9886,
            (
                1.00,
                0.9927,
                0.9927,
                1.00,
                0.9927,
                0.9854,
                0.9927,
                0.9927,
                0.9927,
                0.9794,
                0.9794,
                0.9794,
                0.9794,
            ),
        )

    def test_gc_fraction(self):
        """Tests gc_fraction function."""
        self.assertAlmostEqual(gc_fraction("", "ignore"), 0, places=3)
        self.assertAlmostEqual(gc_fraction("", "weighted"), 0, places=3)
        self.assertAlmostEqual(gc_fraction("", "remove"), 0, places=3)

        seq = "ACGGGCTACCGTATAGGCAAGAGATGATGCCC"
        self.assertAlmostEqual(gc_fraction(seq, "ignore"), 0.5625, places=3)
        self.assertAlmostEqual(gc_fraction(seq, "weighted"), 0.5625, places=3)
        self.assertAlmostEqual(gc_fraction(seq, "remove"), 0.5625, places=3)

        seq = "ACTGSSSS"
        self.assertAlmostEqual(gc_fraction(seq, "ignore"), 0.75, places=3)
        self.assertAlmostEqual(gc_fraction(seq, "weighted"), 0.75, places=3)
        self.assertAlmostEqual(gc_fraction(seq, "remove"), 0.75, places=3)

        # Test RNA sequence

        seq = "GGAUCUUCGGAUCU"
        self.assertAlmostEqual(gc_fraction(seq, "ignore"), 0.5, places=3)
        self.assertAlmostEqual(gc_fraction(seq, "weighted"), 0.5, places=3)
        self.assertAlmostEqual(gc_fraction(seq, "remove"), 0.5, places=3)

        # Test ambiguous nucleotide behaviour

        seq = "CCTGNN"
        self.assertAlmostEqual(gc_fraction(seq, "ignore"), 0.5, places=3)
        self.assertAlmostEqual(gc_fraction(seq, "weighted"), 0.667, places=3)
        self.assertAlmostEqual(gc_fraction(seq, "remove"), 0.75, places=3)

        seq = "GDVV"
        self.assertAlmostEqual(gc_fraction(seq, "ignore"), 0.25, places=3)
        self.assertAlmostEqual(gc_fraction(seq, "weighted"), 0.6667, places=3)
        self.assertAlmostEqual(gc_fraction(seq, "remove"), 1.00, places=3)

        with self.assertRaises(ValueError):
            gc_fraction(seq, "other string")

    def test_GC_skew(self):
        s = "A" * 50
        seq = Seq(s)
        record = SeqRecord(seq)
        self.assertEqual(GC_skew(s)[0], 0)
        self.assertEqual(GC_skew(seq)[0], 0)
        self.assertEqual(GC_skew(record)[0], 0)

    def test_seq1_seq3(self):
        s3 = "MetAlaTyrtrpcysthrLYSLEUILEGlYPrOGlNaSnaLapRoTyRLySSeRHisTrpLysThr"
        s1 = "MAYWCTKLIGPQNAPYKSHWKT"
        self.assertEqual(seq1(s3), s1)
        self.assertEqual(seq3(s1).upper(), s3.upper())
        self.assertEqual(seq1(seq3(s1)), s1)
        self.assertEqual(seq3(seq1(s3)).upper(), s3.upper())

    def test_lcc_simp(self):
        s = "ACGATAGC"
        seq = Seq(s)
        record = SeqRecord(seq)
        self.assertAlmostEqual(lcc_simp(s), 0.9528, places=4)
        self.assertAlmostEqual(lcc_simp(seq), 0.9528, places=4)
        self.assertAlmostEqual(lcc_simp(record), 0.9528, places=4)

    def test_lcc_mult(self):
        s = "ACGATAGC"
        seq = Seq(s)
        record = SeqRecord(seq)
        llc_lst = lcc_mult(s, len(s))
        self.assertEqual(len(llc_lst), 1)
        self.assertAlmostEqual(llc_lst[0], 0.9528, places=4)
        llc_lst = lcc_mult(seq, len(seq))
        self.assertEqual(len(llc_lst), 1)
        self.assertAlmostEqual(llc_lst[0], 0.9528, places=4)
        llc_lst = lcc_mult(record, len(record))
        self.assertEqual(len(llc_lst), 1)
        self.assertAlmostEqual(llc_lst[0], 0.9528, places=4)


class NtSearchTests(unittest.TestCase):
    """Tests for Bio.SeqUtils.nt_search."""

    @staticmethod
    def matching_positions(seq, subseq):
        """Return every position at which subseq matches seq.

        This is deliberately an independent (naive) implementation written
        from the IUPAC ambiguity code definitions, so that its answer can be
        compared against the regular-expression based search in
        ``Bio.SeqUtils.nt_search``.  Note that overlapping matches count, and
        that an ambiguity code in ``seq`` only matches itself (the ambiguity
        codes are expanded on the ``subseq`` side only).
        """
        positions = []
        for start in range(len(seq) - len(subseq) + 1):
            window = seq[start : start + len(subseq)]
            if all(
                base in IUPACData.ambiguous_dna_values[code]
                for base, code in zip(window, subseq)
            ):
                positions.append(start)
        return positions

    def test_agrees_with_naive_search(self):
        """Positions returned agree with a naive independent search."""
        cases = [
            ("AAGATTAGCATCGGATCC", "AT"),
            ("AAAA", "AA"),  # overlapping matches
            ("AAAAA", "AAA"),  # overlapping matches
            ("ACGTACGTACGT", "ACGT"),
            ("GGGTCAGTCAGTCA", "RGT"),  # R = A or G
            ("CATGCATGCATG", "NNN"),  # matches everywhere it fits
            ("ACGTACGT", "Y"),  # Y = C or T
            ("ACGT", "TTTT"),  # subseq longer than any match
            ("ACGT", "GGG"),  # no match at all
            ("TTTTAAAA", "T"),  # run at the very start
            ("AAAATTTT", "T"),  # run at the very end
            ("ACGTNACGT", "N"),  # N in seq is not expanded
        ]
        for seq, subseq in cases:
            with self.subTest(seq=seq, subseq=subseq):
                result = nt_search(seq, subseq)
                self.assertEqual(result[1:], self.matching_positions(seq, subseq))

    def test_returned_pattern_matches_the_hits(self):
        """The regular expression returned first matches every reported hit."""
        for seq, subseq in [("GGGTCAGTCAGTCA", "RGT"), ("ACGTACGTAA", "MRY")]:
            with self.subTest(seq=seq, subseq=subseq):
                result = nt_search(seq, subseq)
                pattern = result[0]
                self.assertIsInstance(pattern, str)
                hits = result[1:]
                self.assertNotEqual(hits, [])
                for position in hits:
                    fragment = seq[position : position + len(subseq)]
                    self.assertIsNotNone(re.fullmatch(pattern, fragment))
                # and it must not match anywhere that was not reported
                for position in range(len(seq) - len(subseq) + 1):
                    if position in hits:
                        continue
                    fragment = seq[position : position + len(subseq)]
                    self.assertIsNone(re.fullmatch(pattern, fragment))

    def test_final_position_is_reported(self):
        """A match ending at the last base is reported (off-by-one guard)."""
        self.assertEqual(nt_search("CCCCAT", "AT")[1:], [4])
        self.assertEqual(nt_search("AT", "AT")[1:], [0])

    def test_unknown_letter_in_subseq(self):
        """An illegal letter in the subsequence raises KeyError."""
        with self.assertRaises(KeyError):
            nt_search("ACGT", "AZ")

    def test_empty_subseq_is_rejected(self):
        """An empty subsequence raises ValueError instead of looping forever."""
        for seq in ("AAGATTAGCATCGGATCC", "", Seq("ACGT")):
            for subseq in ("", Seq(""), MutableSeq("")):
                with self.subTest(seq=str(seq), subseq=type(subseq).__name__):
                    with self.assertRaises(ValueError) as context:
                        nt_search(seq, subseq)
                    self.assertEqual(str(context.exception), "subseq must not be empty")

    def test_accepts_sequence_objects(self):
        """Seq, MutableSeq and SeqRecord are searched like the equivalent string."""
        text = "AAGATTAGCATCGGATCC"
        # AT starts at index 3, at index 9 (A at 9, T at 10) and at index 14
        expected = ["AT", 3, 9, 14]
        self.assertEqual(nt_search(text, "AT"), expected)
        for sequence in (Seq(text), MutableSeq(text), SeqRecord(Seq(text), id="x")):
            with self.subTest(sequence=type(sequence).__name__):
                self.assertEqual(nt_search(sequence, "AT"), expected)
        for subseq in (Seq("AT"), MutableSeq("AT")):
            with self.subTest(subseq=type(subseq).__name__):
                self.assertEqual(nt_search(text, subseq), expected)


class GC123Tests(unittest.TestCase):
    """Tests for Bio.SeqUtils.GC123."""

    def test_by_codon_position(self):
        """G+C is reported overall and for each codon position."""
        # ATG GCC: overall 4 of 6 are G or C; first positions A, G -> 1 of 2;
        # second positions T, C -> 1 of 2; third positions G, C -> 2 of 2.
        overall, first, second, third = GC123("ATGGCC")
        self.assertAlmostEqual(overall, 400 / 6)
        self.assertEqual((first, second, third), (50.0, 50.0, 100.0))

    def test_ambiguous_nucleotides_are_not_counted(self):
        """Only A, C, G and T count towards the totals."""
        # A C T G T N: five countable bases of which C and G are G+C, so 40%;
        # first positions A, G -> 50%, second C, T -> 50%, third T only -> 0%.
        self.assertEqual(GC123("ACTGTN"), (40.0, 50.0, 50.0, 0.0))

    def test_codon_position_without_any_nucleotide(self):
        """A codon position with nothing countable scores zero, not an error."""
        # A C N: the third position has no countable base at all.
        self.assertEqual(GC123("ACN"), (50.0, 0.0, 100.0, 0.0))

    def test_incomplete_final_codon(self):
        """A trailing partial codon is padded rather than dropped."""
        # ACT G: the G is a first codon position.
        self.assertEqual(GC123("ACTG"), (50.0, 50.0, 100.0, 0.0))

    def test_case_insensitive(self):
        """Lower case sequences give the same answer."""
        self.assertEqual(GC123("atggcc"), GC123("ATGGCC"))


class Seq1Seq3EdgeCaseTests(unittest.TestCase):
    """Edge cases of Bio.SeqUtils.seq1 and Bio.SeqUtils.seq3."""

    def test_seq3_does_not_modify_iupac_tables(self):
        """seq3 must not mutate IUPACData when given a custom map."""
        before = dict(IUPACData.protein_letters_1to3_extended)
        self.assertEqual(seq3("A*", custom_map={"*": "***"}), "Ala***")
        self.assertEqual(IUPACData.protein_letters_1to3_extended, before)

    def test_seq1_does_not_modify_iupac_tables(self):
        """seq1 must not mutate IUPACData when given a custom map."""
        before = dict(IUPACData.protein_letters_3to1_extended)
        self.assertEqual(seq1("AlaTer", custom_map={"Ter": "#"}), "A#")
        self.assertEqual(IUPACData.protein_letters_3to1_extended, before)

    def test_seq1_ignores_trailing_partial_codon(self):
        """seq1 reads whole three-letter blocks and drops any remainder."""
        self.assertEqual(seq1("MetAla"), "MA")
        self.assertEqual(seq1("MetAlaX"), "MA")
        self.assertEqual(seq1("MetAlaXX"), "MA")
        self.assertEqual(seq1("MetAlaGly"), "MAG")

    def test_seq1_custom_map_replaces_default_terminator(self):
        """A custom map overrides the default {'Ter': '*'} mapping."""
        self.assertEqual(seq1("MetAlaTer"), "MA*")
        self.assertEqual(seq1("MetAlaTer", custom_map={"Ter": "X"}), "MAX")

    def test_seq1_custom_map_is_case_insensitive(self):
        """Custom map keys are matched case-insensitively, like the defaults."""
        self.assertEqual(seq1("metalater", custom_map={"ter": "+"}), "MA+")
        self.assertEqual(seq1("METALATER", custom_map={"Ter": "+"}), "MA+")

    def test_seq3_undef_code(self):
        """Unknown one-letter codes become undef_code."""
        self.assertEqual(seq3("A-B"), "AlaXaaAsx")
        self.assertEqual(seq3("A-B", undef_code="???"), "Ala???Asx")

    def test_seq1_undef_code(self):
        """Unknown three-letter codes become undef_code."""
        self.assertEqual(seq1("AlaZzzGly"), "AXG")
        self.assertEqual(seq1("AlaZzzGly", undef_code="?"), "A?G")


class MolecularWeightTests(unittest.TestCase):
    """Tests for Bio.SeqUtils.molecular_weight."""

    # Average and monoisotopic masses of water.  The monoisotopic value is
    # 2 x 1.0078250319 (1H) + 15.9949146221 (16O).
    water = 18.0153
    monoisotopic_water = 2 * 1.0078250319 + 15.9949146221

    def test_dna_from_published_nucleotide_masses(self):
        """DNA mass equals sum of dNMP masses minus one water per bond.

        The average masses of the free-acid 2'-deoxyribonucleoside
        5'-monophosphates are dAMP 331.2218, dCMP 307.1971, dGMP 347.2212
        and dTMP 322.2085 g/mol; each phosphodiester bond releases one
        molecule of water (18.0153 g/mol).
        """
        dAMP, dCMP, dGMP, dTMP = 331.2218, 307.1971, 347.2212, 322.2085
        expected = dAMP + dGMP + dCMP - 2 * self.water
        self.assertAlmostEqual(molecular_weight("AGC", "DNA"), expected, places=4)
        expected = dAMP + dCMP + dGMP + dTMP - 3 * self.water
        self.assertAlmostEqual(molecular_weight("ACGT", "DNA"), expected, places=4)

    def test_protein_from_published_amino_acid_masses(self):
        """Peptide mass equals sum of amino acid masses minus one water per bond.

        Average masses of the free amino acids: Ala 89.0932, Gly 75.0666,
        Cys 121.1582 g/mol.
        """
        ala, gly, cys = 89.0932, 75.0666, 121.1582
        expected = ala + gly + cys - 2 * self.water
        self.assertAlmostEqual(molecular_weight("AGC", "protein"), expected, places=4)

    def test_concatenation_loses_exactly_one_water(self):
        """MW(a + b) == MW(a) + MW(b) - water, for every sequence type.

        This pins the ``(len(seq) - 1) * water`` term: any other multiplier
        makes the identity fail.
        """
        cases = [
            ("DNA", "ACGTTGCA", "TTGACCAG", False),
            ("RNA", "ACGUUGCA", "UUGACCAG", False),
            ("protein", "MKWVTFISL", "LLFSSAYSR", False),
            ("DNA", "ACGTTGCA", "TTGACCAG", True),
            ("protein", "MKWVTFISL", "LLFSSAYSR", True),
        ]
        for seq_type, first, second, monoisotopic in cases:
            with self.subTest(seq_type=seq_type, monoisotopic=monoisotopic):
                water = self.monoisotopic_water if monoisotopic else self.water
                joined = molecular_weight(
                    first + second, seq_type, monoisotopic=monoisotopic
                )
                separate = molecular_weight(
                    first, seq_type, monoisotopic=monoisotopic
                ) + molecular_weight(second, seq_type, monoisotopic=monoisotopic)
                self.assertAlmostEqual(joined, separate - water, places=4)

    def test_single_residue_has_no_water_loss(self):
        """A one residue sequence weighs exactly one residue."""
        for letter in "ACGT":
            with self.subTest(letter=letter):
                self.assertAlmostEqual(
                    molecular_weight(letter + letter, "DNA"),
                    2 * molecular_weight(letter, "DNA") - self.water,
                    places=4,
                )

    def test_monoisotopic_water_constant(self):
        """The monoisotopic water constant matches the atomic masses."""
        difference = 2 * molecular_weight("A", "DNA", monoisotopic=True) - (
            molecular_weight("AA", "DNA", monoisotopic=True)
        )
        self.assertAlmostEqual(difference, self.monoisotopic_water, places=5)

    def test_monoisotopic_is_lighter_than_average(self):
        """Monoisotopic masses are below the average masses."""
        for seq_type, seq in [
            ("DNA", "ACGTTGCA"),
            ("RNA", "ACGUUGCA"),
            ("protein", "MKWVTFISL"),
        ]:
            with self.subTest(seq_type=seq_type):
                self.assertLess(
                    molecular_weight(seq, seq_type, monoisotopic=True),
                    molecular_weight(seq, seq_type, monoisotopic=False),
                )

    def test_circular_loses_one_more_water(self):
        """A circular molecule weighs one water less than the linear one."""
        for seq_type, seq in [
            ("DNA", "ACGTTGCA"),
            ("RNA", "ACGUUGCA"),
            ("protein", "MKWVTFISL"),
        ]:
            with self.subTest(seq_type=seq_type):
                self.assertAlmostEqual(
                    molecular_weight(seq, seq_type, circular=True),
                    molecular_weight(seq, seq_type) - self.water,
                    places=4,
                )

    def test_double_stranded_is_sum_of_both_strands(self):
        """Double stranded mass equals the mass of both single strands."""
        seq = "ACGTTGCAA"
        self.assertAlmostEqual(
            molecular_weight(seq, "DNA", double_stranded=True),
            molecular_weight(seq, "DNA") + molecular_weight(complement(seq), "DNA"),
            places=4,
        )
        rna = "ACGUUGCAA"
        self.assertAlmostEqual(
            molecular_weight(rna, "RNA", double_stranded=True),
            molecular_weight(rna, "RNA") + molecular_weight(complement_rna(rna), "RNA"),
            places=4,
        )

    def test_double_stranded_is_strand_symmetric(self):
        """A duplex weighs the same whichever strand is given."""
        seq = "ACGTTGCAAGGTC"
        self.assertAlmostEqual(
            molecular_weight(seq, "DNA", double_stranded=True),
            molecular_weight(reverse_complement(seq), "DNA", double_stranded=True),
            places=4,
        )

    def test_double_stranded_circular(self):
        """A circular duplex loses one water per strand."""
        seq = "ACGTTGCAAGGTC"
        self.assertAlmostEqual(
            molecular_weight(seq, "DNA", double_stranded=True, circular=True),
            molecular_weight(seq, "DNA", double_stranded=True) - 2 * self.water,
            places=4,
        )

    def test_accepts_seq_and_seqrecord_and_whitespace(self):
        """Seq, SeqRecord, whitespace and lower case all give the same answer."""
        expected = molecular_weight("ACGTTGCA", "DNA")
        for value in (
            Seq("ACGTTGCA"),
            MutableSeq("ACGTTGCA"),
            SeqRecord(Seq("ACGTTGCA")),
            "acgttgca",
            "ACGT TGCA",
            "ACGT\nTGCA",
        ):
            with self.subTest(value=repr(value)[:40]):
                self.assertAlmostEqual(
                    molecular_weight(value, "DNA"), expected, places=4
                )

    def test_bad_sequence_type(self):
        """An unknown sequence type is rejected."""
        with self.assertRaises(ValueError) as context:
            molecular_weight("ACGT", "peptide")
        self.assertIn("peptide", str(context.exception))

    def test_ambiguous_letter_rejected(self):
        """Ambiguous or unknown letters are rejected."""
        for seq_type, seq in [("DNA", "ACNT"), ("RNA", "ACNU"), ("protein", "MKZ")]:
            with self.subTest(seq_type=seq_type):
                with self.assertRaises(ValueError) as context:
                    molecular_weight(seq, seq_type)
                self.assertIn(seq_type, str(context.exception))

    def test_protein_cannot_be_double_stranded(self):
        """Proteins have no complement."""
        with self.assertRaises(ValueError):
            molecular_weight("MKWV", "protein", double_stranded=True)


class SixFrameTranslationsTests(unittest.TestCase):
    """Tests for Bio.SeqUtils.six_frame_translations."""

    @staticmethod
    def frames(seq):
        """Return the six translations, computed independently.

        Frame ``+n`` starts at offset ``n - 1`` of ``seq``; frame ``-n``
        starts at offset ``n - 1`` of the reverse complement.  The reverse
        frames are reported 3' to 5' with respect to ``seq``, i.e. reversed,
        which is the convention ``six_frame_translations`` uses.
        """
        result = {}
        antiparallel = reverse_complement(seq)
        for offset in range(3):
            forward = [seq[i : i + 3] for i in range(offset, len(seq) - 2, 3)]
            reverse = [
                antiparallel[i : i + 3] for i in range(offset, len(antiparallel) - 2, 3)
            ]
            result[offset + 1] = "".join(translate(codon) for codon in forward)
            result[-(offset + 1)] = "".join(
                translate(codon) for codon in reversed(reverse)
            )
        return result

    def test_short_dna_sequence(self):
        """Header and frames of a short (<= 20 nt) DNA sequence."""
        seq = "ATGGCCATTGTA"
        lines = six_frame_translations(seq).split("\n")
        # 3 A, 4 T, 3 G and 2 C; GC content 5/12 = 41.67%
        self.assertEqual(lines[0], "GC_Frame: a:3 t:4 g:3 c:2")
        # Sequences of 20 nt or less are shown in full, not abbreviated
        self.assertEqual(lines[1], "Sequence: atggccattgta, 12 nt, 41.67 %GC")
        self.assertEqual(lines[4], "1/1")
        self.assertEqual(lines[8].split()[0], seq.lower())
        self.assertEqual(lines[9], complement(seq).lower())

    def test_long_sequence_header_is_abbreviated(self):
        """Sequences longer than 20 nt are abbreviated in the header."""
        seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
        lines = six_frame_translations(seq).split("\n")
        gc = 100 * (seq.count("G") + seq.count("C")) / len(seq)
        self.assertEqual(
            lines[1],
            "Sequence: %s ... %s, %d nt, %0.2f %%GC"
            % (seq[:10].lower(), seq[-10:].lower(), len(seq), gc),
        )

    def test_forward_frames_are_correctly_placed(self):
        """Each forward frame residue sits above the codon it comes from."""
        seq = "ATGGCCATTGTAATGGGCCGCTGA"
        lines = six_frame_translations(seq).split("\n")
        # rows 5, 6 and 7 hold frames +3, +2 and +1 respectively
        for row, frame in ((5, 3), (6, 2), (7, 1)):
            with self.subTest(frame=frame):
                text = lines[row]
                for column, residue in enumerate(text):
                    if residue == " ":
                        continue
                    codon = seq[column : column + 3]
                    self.assertEqual(
                        residue,
                        translate(codon),
                        f"frame {frame}, column {column}, codon {codon}",
                    )

    def test_all_six_frames_are_present(self):
        """The output holds all six independently computed translations.

        The three forward frames are printed above the sequence and the
        three reverse ones below it.  Only the content is checked here, not
        which reverse frame is printed on which row.
        """
        for seq in (
            "ATGGCCATTGTAATGGGCCGCTGA",  # length a multiple of three
            "ATGGCCATTGTAATGGGCCGCTGAA",  # and the two other cases
            "ATGGCCATTGTAATGGGCCGCTGAAC",
        ):
            with self.subTest(seq=seq):
                lines = six_frame_translations(seq).split("\n")
                expected = self.frames(seq)
                forward = [lines[row].replace(" ", "") for row in (5, 6, 7)]
                reverse = [lines[row].replace(" ", "") for row in (10, 11, 12)]
                self.assertEqual(forward, [expected[3], expected[2], expected[1]])
                self.assertEqual(
                    sorted(reverse),
                    sorted([expected[-1], expected[-2], expected[-3]]),
                )

    def test_reverse_frames_are_correctly_placed(self):
        """Each reverse frame residue sits below the codon it comes from.

        A residue printed at column ``c`` of a reverse-strand row is the
        translation of the reverse complement of ``seq[c:c + 3]``.  Frame
        ``-(n + 1)`` is translated from offset ``n`` of the reverse
        complement, so its leftmost residue lands at column
        ``(len(seq) - n) % 3`` of the top strand; the layout therefore
        depends on ``len(seq) % 3``, and all three cases are covered here.
        """
        for seq in (
            "ATGGCCATTGTAATGGGCCGCTGA",  # 24 nt, 24 % 3 == 0
            "ATGGCCATTGTAATGGGCCGCTGAA",  # 25 nt, 25 % 3 == 1
            "ATGGCCATTGTAATGGGCCGCTGAAC",  # 26 nt, 26 % 3 == 2
        ):
            with self.subTest(seq=seq):
                lines = six_frame_translations(seq).split("\n")
                # rows 10, 11 and 12 hold the three reverse frames
                for row in (10, 11, 12):
                    text = lines[row]
                    self.assertNotEqual(text.strip(), "")
                    for column, residue in enumerate(text):
                        if residue == " ":
                            continue
                        codon = seq[column : column + 3]
                        self.assertEqual(
                            residue,
                            translate(reverse_complement(codon)),
                            f"row {row}, column {column}, codon {codon}",
                        )

    def test_reverse_frame_columns_of_a_hand_worked_example(self):
        """The three reverse rows of a worked example, derived by hand.

        ``ATGAAATTTGGG`` is 12 nt, so ``len(seq) % 3 == 0`` and frame
        ``-(n + 1)`` starts at column ``(12 - n) % 3``: frame -1 at column 0,
        frame -2 at column 2 and frame -3 at column 1.  Reverse complementing
        the top-strand codons at those columns gives, for frame -1,
        ATG AAA TTT GGG -> CAT TTT AAA CCC -> H F K P; for frame -3,
        TGA AAT TTG -> TCA ATT CAA -> S I Q; and for frame -2,
        GAA ATT TGG -> TTC AAT CCA -> F N P.
        """
        lines = six_frame_translations("ATGAAATTTGGG").split("\n")
        self.assertEqual(lines[10], "H  F  K  P")
        self.assertEqual(lines[11], " S  I  Q")
        self.assertEqual(lines[12], "  F  N  P")

    def test_rna_input(self):
        """An RNA sequence is reverse complemented as RNA."""
        rna = "AUGGCCAUUGUA"
        lines = six_frame_translations(rna).split("\n")
        self.assertEqual(lines[9], complement_rna(rna).lower())
        self.assertEqual(lines[7], "  ".join(translate(rna)))


class CodonAdaptationIndexErrorTests(unittest.TestCase):
    """Error handling and small hand-checked cases for CodonAdaptationIndex."""

    def test_relative_adaptiveness_from_counts(self):
        """The w values follow Sharp & Li's definition.

        Sharp & Li (1987) define the relative adaptiveness of a codon as its
        observed frequency divided by that of the most frequently used
        synonymous codon; codons that were never seen are given a count of
        0.5.  The reference gene here is ATG AAA AAG AAA TAA, so the counts
        are ATG 1, AAA 2, AAG 1, TAA 1 and 0.5 everywhere else.
        """
        index = CodonAdaptationIndex(["ATGAAAAAGAAATAA"])
        self.assertAlmostEqual(index["AAA"], 1.0)  # 2 / 2
        self.assertAlmostEqual(index["AAG"], 0.5)  # 1 / 2
        self.assertAlmostEqual(index["ATG"], 1.0)  # only Met codon
        self.assertAlmostEqual(index["TAA"], 1.0)  # 1 / 1
        self.assertAlmostEqual(index["TAG"], 0.5)  # 0.5 / 1
        self.assertAlmostEqual(index["TGA"], 0.5)  # 0.5 / 1
        # Phe was never seen, so both its codons get 0.5 / 0.5 == 1
        self.assertAlmostEqual(index["TTT"], 1.0)
        self.assertAlmostEqual(index["TTC"], 1.0)

    def test_calculate_is_geometric_mean(self):
        """CAI is the geometric mean of the w values of the counted codons."""
        index = CodonAdaptationIndex(["ATGAAAAAGAAATAA"])
        # AAG (w = 0.5) and AAA (w = 1.0): sqrt(0.5 * 1.0)
        self.assertAlmostEqual(index.calculate("AAGAAA"), 0.5**0.5, places=10)
        # ATG and TGG are excluded, so only AAG is counted here
        self.assertAlmostEqual(index.calculate("ATGAAG"), 0.5, places=10)
        self.assertAlmostEqual(index.calculate("AAGTGG"), 0.5, places=10)

    def test_illegal_codon_without_gene_name(self):
        """A bad codon in a bare sequence is reported without a gene name."""
        for sequence in ("ATGAANAAA", Seq("ATGAANAAA"), MutableSeq("ATGAANAAA")):
            with self.subTest(sequence=type(sequence).__name__):
                with self.assertRaises(ValueError) as context:
                    CodonAdaptationIndex([sequence])
                self.assertEqual(str(context.exception), "illegal codon 'AAN'")

    def test_illegal_codon_with_gene_name(self):
        """A bad codon in a SeqRecord is reported with the record id."""
        record = SeqRecord(Seq("ATGAANAAA"), id="test_gene")
        with self.assertRaises(ValueError) as context:
            CodonAdaptationIndex([record])
        self.assertEqual(
            str(context.exception), "illegal codon 'AAN' in gene test_gene"
        )

    def test_mutable_seq_sequences_are_accepted(self):
        """MutableSeq input works, as the docstrings say it should.

        Slicing a MutableSeq gives a MutableSeq, which is unhashable and so
        cannot be used to look a codon up in the count dictionary.
        """
        reference = "ATGAAAAAGAAATAA"
        expected = CodonAdaptationIndex([reference])
        for sequence in (
            MutableSeq(reference),
            SeqRecord(MutableSeq(reference), id="test_gene"),
        ):
            with self.subTest(sequence=type(sequence).__name__):
                self.assertEqual(dict(CodonAdaptationIndex([sequence])), dict(expected))
        # and calculate() accepts one too: AAG (w = 0.5) and AAA (w = 1.0)
        self.assertAlmostEqual(
            expected.calculate(MutableSeq("AAGAAA")), 0.5**0.5, places=10
        )
        self.assertAlmostEqual(
            expected.calculate(SeqRecord(MutableSeq("AAGAAA"))), 0.5**0.5, places=10
        )

    def test_lower_case_sequences_are_accepted(self):
        """Sequences are upper cased before counting."""
        upper = CodonAdaptationIndex(["ATGAAAAAGAAATAA"])
        lower = CodonAdaptationIndex(["atgaaaaagaaataa"])
        self.assertEqual(dict(upper), dict(lower))

    def test_calculate_skips_missing_stop_codons(self):
        """A stop codon missing from the index is skipped, not an error."""
        index = CodonAdaptationIndex(["ATGAAAAAGAAATAA"])
        for codon in ("TAA", "TAG", "TGA"):
            del index[codon]
        self.assertAlmostEqual(index.calculate("AAGTAA"), 0.5, places=10)
        self.assertAlmostEqual(index.calculate("AAGTAG"), 0.5, places=10)
        self.assertAlmostEqual(index.calculate("AAGTGA"), 0.5, places=10)

    def test_calculate_rejects_missing_sense_codon(self):
        """A sense codon missing from the index is an error."""
        index = CodonAdaptationIndex(["ATGAAAAAGAAATAA"])
        del index["AAA"]
        with self.assertRaises(TypeError) as context:
            index.calculate("AAGAAA")
        self.assertEqual(str(context.exception), "illegal codon in sequence: AAA")

    def test_optimize_rejects_ties_when_strict(self):
        """Two equally preferred codons are an error when strict."""
        # Without reference sequences every codon has w == 1.0, so every
        # amino acid with more than one codon is ambiguous.
        index = CodonAdaptationIndex([])
        with self.assertRaises(ValueError) as context:
            index.optimize("AAAAAG")
        self.assertIn("equally preferred", str(context.exception))

    def test_optimize_rejects_unknown_sequence_type(self):
        """An unknown seq_type is rejected by optimize."""
        index = CodonAdaptationIndex([])
        with self.assertRaises(ValueError) as context:
            index.optimize("AAAAAG", seq_type="peptide", strict=False)
        self.assertIn("peptide", str(context.exception))

    def test_optimize_preserves_the_protein(self):
        """Optimising does not change the encoded protein."""
        index = CodonAdaptationIndex(["ATGAAAAAGAAATAA"])
        for seq_type, sequence in (
            ("DNA", "ATGAAGAAA"),
            ("RNA", "AUGAAGAAA"),
            ("protein", "MKK"),
        ):
            with self.subTest(seq_type=seq_type):
                optimized = index.optimize(sequence, seq_type, strict=False)
                self.assertEqual(translate(optimized), "MKK")
                self.assertEqual(index.calculate(optimized), 1.0)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
