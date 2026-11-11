# Copyright (C) 2013 by Zheng Ruan (zruan1991@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the Bio.codonalign modules."""

import math
import os
import subprocess
import sys
import tempfile
import unittest
import warnings

import support

from Bio import AlignIO
from Bio import BiopythonWarning
from Bio import codonalign
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

CODONALIGN = support.DATA / "codonalign"
TEST_ALIGN_FILE1 = [(CODONALIGN / "nucl1.fa", CODONALIGN / "pro1.aln"), "parse"]
TEST_ALIGN_FILE2 = [(CODONALIGN / "nucl2.fa", CODONALIGN / "pro2.aln"), "parse"]
TEST_ALIGN_FILE3 = [(CODONALIGN / "nucl3.fa", CODONALIGN / "pro3.aln"), "index"]
TEST_ALIGN_FILE4 = [(CODONALIGN / "nucl4.fa", CODONALIGN / "pro4.aln"), "index"]
TEST_ALIGN_FILE5 = [(CODONALIGN / "nucl5.fa", CODONALIGN / "pro5.aln"), "parse"]
TEST_ALIGN_FILE6 = [
    (CODONALIGN / "egfr_nucl.fa", CODONALIGN / "egfr_pro.aln", CODONALIGN / "egfr_id"),
    "id",
]
TEST_ALIGN_FILE7 = [(CODONALIGN / "drosophila.fasta", CODONALIGN / "adh.aln"), "index"]

temp_dir = tempfile.mkdtemp()


class TestImport(unittest.TestCase):
    def test_no_experimental_warning(self):
        """Importing Bio.codonalign does not warn that it is experimental.

        Run in a subprocess because this process has already imported it.
        """
        code = (
            "import warnings\n"
            "from Bio import BiopythonExperimentalWarning\n"
            "warnings.simplefilter('error', BiopythonExperimentalWarning)\n"
            "import Bio.codonalign\n"
        )
        result = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
            cwd=os.path.dirname(os.path.abspath(__file__)),
        )
        self.assertEqual(result.returncode, 0, result.stderr)


class TestCodonSeq(unittest.TestCase):
    def test_seq(self):
        codonseq1 = codonalign.CodonSeq("AAATTT---TTTGGACCC", rf_table=[0, 3, 6, 9, 12])
        self.assertEqual(len(codonseq1), 18)
        self.assertEqual(codonseq1.get_codon_num(), 5)
        self.assertEqual(str(codonseq1.get_codon(0)), "AAA")
        self.assertEqual(str(codonseq1.get_codon(-1)), "CCC")
        self.assertEqual(str(codonseq1.get_codon(slice(1, 3))), "TTT---")
        self.assertEqual(
            str(codonseq1.get_codon(slice(None, None, -1))), "CCCGGATTT---TTTAAA"
        )

        self.assertRaises(ValueError, codonalign.CodonSeq, "AAA-T")
        self.assertIsInstance(codonseq1.toSeq(), Seq)


class TestCodonAlignment(unittest.TestCase):
    def setUp(self):
        codonseq1 = codonalign.CodonSeq("AAATTT---TTTGGACCC")
        codonseq2 = codonalign.CodonSeq("AAGTTT---TTTGGGCCC")
        codonseq3 = codonalign.CodonSeq("AAGTAT---TTTGGACCC")
        codonseq4 = codonalign.CodonSeq("AACTTT---TTTGGACGC")

        self.seqrec = [
            SeqRecord(codonseq1, id="alpha"),
            SeqRecord(codonseq2, id="beta"),
            SeqRecord(codonseq3, id="gamma"),
            SeqRecord(codonseq4, id="delta"),
        ]

    def test_align(self):
        codonAlign = codonalign.CodonAlignment(self.seqrec)
        self.assertEqual(codonAlign.get_aln_length(), 6)
        self.assertIsInstance(codonAlign.toMultipleSeqAlignment(), MultipleSeqAlignment)


class TestAddition(unittest.TestCase):
    def setUp(self):
        self.seq1 = SeqRecord(Seq("ATGTCTCGT"), id="pro1")
        self.seq2 = SeqRecord(Seq("ATGCGT"), id="pro2")
        self.pro1 = SeqRecord(Seq("MSR"), id="pro1")
        self.pro2 = SeqRecord(Seq("M-R"), id="pro2")
        self.aln = MultipleSeqAlignment([self.pro1, self.pro2])
        self.codon_aln = codonalign.build(self.aln, [self.seq1, self.seq2])

        tail1 = SeqRecord(Seq("AAA"), id="pro1")
        tail2 = SeqRecord(Seq("AAA"), id="pro2")
        self.multi_aln = MultipleSeqAlignment([tail1, tail2])

    def test_addition_MultipleSeqAlignment(self):
        """Check addition of CodonAlignment and MultipleSeqAlignment."""
        new_aln1 = self.codon_aln + self.multi_aln

        self.assertIsInstance(new_aln1, MultipleSeqAlignment)
        for x in range(len(self.codon_aln)):
            self.assertEqual(
                new_aln1[x].seq, self.codon_aln[x].seq + self.multi_aln[x].seq
            )

        new_aln2 = self.multi_aln + self.codon_aln

        self.assertIsInstance(new_aln2, MultipleSeqAlignment)
        for x in range(len(self.codon_aln)):
            self.assertEqual(
                new_aln2[x].seq, self.multi_aln[x].seq + self.codon_aln[x].seq
            )

    def test_addition_CodonAlignment(self):
        """Check addition of CodonAlignment and CodonAlignment."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=BiopythonWarning)
            new_aln = self.codon_aln + self.codon_aln

        self.assertIsInstance(new_aln, codonalign.CodonAlignment)
        for x in range(len(self.codon_aln)):
            self.assertEqual(
                new_aln[x].seq, self.codon_aln[x].seq + self.codon_aln[x].seq
            )

    def test_ValueError(self):
        """Check that ValueError is thrown for Alignments of different lengths."""
        # original len(self.aln) = 2 , len(aln) = 3
        aln = MultipleSeqAlignment(
            [self.pro1, self.pro2, SeqRecord(Seq("M--"), id="pro3")]
        )
        triple_codon = codonalign.build(
            aln, [self.seq1, self.seq2, SeqRecord(Seq("ATG"), id="pro3")]
        )
        with self.assertRaises(ValueError):
            triple_codon + self.multi_aln
        with self.assertRaises(ValueError):
            triple_codon + self.codon_aln

    def test_TypeError(self):
        """Check that TypeError is thrown for non CodonAlignment/MultipleSequenceAlignment objects."""
        for obj in [0, "string", ["str1", "str2"], Seq("ATGTCTCGT")]:
            with self.assertRaises(TypeError):
                self.codon_aln + obj


class TestBuildAndIO(unittest.TestCase):
    def setUp(self):
        self.aln_file = [
            TEST_ALIGN_FILE1,
            TEST_ALIGN_FILE2,
            TEST_ALIGN_FILE3,
            TEST_ALIGN_FILE4,
            TEST_ALIGN_FILE5,
            TEST_ALIGN_FILE6,
        ]
        alns = []
        for i in self.aln_file:
            if i[1] == "parse":
                nucl = SeqIO.parse(i[0][0], "fasta")
                prot = AlignIO.read(i[0][1], "clustal")
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    caln = codonalign.build(prot, nucl)
            elif i[1] == "index":
                nucl = SeqIO.index(i[0][0], "fasta")
                prot = AlignIO.read(i[0][1], "clustal")
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    caln = codonalign.build(prot, nucl, max_score=20)
                nucl.close()  # Close the indexed FASTA file
            elif i[1] == "id":
                nucl = SeqIO.parse(i[0][0], "fasta")
                prot = AlignIO.read(i[0][1], "clustal")
                with open(i[0][2]) as handle:
                    id = {i.split()[0]: i.split()[1] for i in handle}
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    caln = codonalign.build(prot, nucl, corr_dict=id)
            alns.append(caln)
        self.alns = alns

    def test_IO(self):
        self.assertEqual(len(self.alns), 6)
        for n, i in enumerate(self.alns):
            aln = i.toMultipleSeqAlignment()
            AlignIO.write(aln, temp_dir + "/aln" + str(n) + ".clw", "clustal")


class Test_build(unittest.TestCase):
    def setUp(self):
        # Test set 1
        seq1 = SeqRecord(
            Seq(
                "TCAGGGACTGCGAGAACCAAGCTACTGCTGCTGCTGGCTGCGCTCTGCGCCGCAGGTGGGGCGCTGGAG"
            ),
            id="pro1",
        )
        seq2 = SeqRecord(
            Seq("TCAGGGACTTCGAGAACCAAGCGCTCCTGCTGCTGGCTGCGCTCGGCGCCGCAGGTGGAGCACTGGAG"),
            id="pro2",
        )
        pro1 = SeqRecord(Seq("SGTARTKLLLLLAALCAAGGALE"), id="pro1")
        pro2 = SeqRecord(Seq("SGTSRTKRLLLLAALGAAGGALE"), id="pro2")
        aln1 = MultipleSeqAlignment([pro1, pro2])
        self.aln1 = aln1
        self.seqlist1 = [seq1, seq2]
        # Test set 2
        #                      M  K  K  H  E L(F)L  C  Q  G  T  S  N  K  L  T  Q(L)L  G  T  F  E  D  H  F  L  S  L  Q  R  M  F  N  N  C  E  V  V
        seq3 = SeqRecord(
            Seq(
                "ATGAAAAAGCACGAGTTACTTTGCCAAGGGACAAGTAACAAGCTCACCCAGTTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAACTGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCTTTTCTCAAGACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC"
            ),
            id="pro1",
        )
        # seq4 =SeqRecord(Seq('ATGAAAAAGCACGAGTT CTTTGCCAAGGGACAAGTAACAAGCTCACCCAGTTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAA TGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCTTTTCTCAAGACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC'), id='pro2')
        seq4 = SeqRecord(
            Seq(
                "ATGAAAAAGCACGAGTTCTTTGCCAAGGGACAAGTAACAAGCTCACCCAGTTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAATGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCTTTTCTCAAGACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC"
            ),
            id="pro2",
        )
        # seq5 =SeqRecord(Seq('ATGAAAAAGCACGAGTT CTTTGCCAAGGGACAAGTAACAAGCTCACCC  TTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAACTGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCTTTTCTCAAGACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC'), id='pro3')
        seq5 = SeqRecord(
            Seq(
                "ATGAAAAAGCACGAGTTACTTTGCCAAGGGACAAGTAACAAGCTCACCCTTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAACTGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCTTTTCTCAAGACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC"
            ),
            id="pro3",
        )
        pro3 = SeqRecord(
            Seq(
                "MKKHELLCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYMQSSYNLSFLKTIQEVAGYVLIAL"
            ),
            id="pro1",
        )
        pro4 = SeqRecord(
            Seq(
                "MKKHEFLCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYMQSSYNLSFLKTIQEVAGYVLIAL"
            ),
            id="pro2",
        )
        pro5 = SeqRecord(
            Seq(
                "MKKHELLCQGTSNKLTLLGTFEDHFLSLQRMFNNCEVVLGNLEITYMQSSYNLSFLKTIQEVAGYVLIAL"
            ),
            id="pro3",
        )
        aln2 = MultipleSeqAlignment([pro3, pro4, pro5])
        self.aln2 = aln2
        self.seqlist2 = [seq3, seq4, seq5]

        # Test set 3
        # use Yeast mitochondrial codon table
        seq6 = SeqRecord(
            Seq(
                "ATGGCAAGGGACCACCCAGTTGGGCACTGATATGATCGGGTGTATTTGCAGAGTAGTAACCTTTCTTTTCTCAAGACCATCCAG"
            ),
            id="pro6",
        )
        seq7 = SeqRecord(
            Seq(
                "ATGGCAAGGCACCATCCAGTTGAGCACTGATATGATCGGGTGTATTTGCAGAGTAGTAACGTGTCTCTGCTCAAGACCATCCAG"
            ),
            id="pro7",
        )
        seq8 = SeqRecord(
            Seq(
                "ATGGCAGGGGACCACCCAGTTGGGCACTGATATGATCGTGTGTATCTGCAGAGTAGTAACCACTCTTTTCTCATGACCATCCAG"
            ),
            id="pro8",
        )
        pro6 = SeqRecord(Seq("MARDHPVGHWYDRVYLQSSNTSFTKTIQ"), id="pro6")
        pro7 = SeqRecord(Seq("MARHHPVEHWYDRVYLQSSNVSTTKTIQ"), id="pro7")
        pro8 = SeqRecord(Seq("MAGDHPVGHWYDRVYTQSSNHSFTMTIQ"), id="pro8")
        aln3 = MultipleSeqAlignment([pro6, pro7, pro8])
        self.aln3 = aln3
        self.seqlist3 = [seq6, seq7, seq8]
        self.codontable3 = CodonTable.unambiguous_dna_by_id[3]

    def test_build(self):
        codon_aln1 = codonalign.build(self.aln1, self.seqlist1)
        codon_aln2 = codonalign.build(self.aln2, self.seqlist2)
        codon_aln3 = codonalign.build(
            self.aln3, self.seqlist3, codon_table=self.codontable3
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=BiopythonWarning)
            codon_aln4 = codonalign.build(
                self.aln1, self.seqlist1, complete_protein=True
            )

    def test_build_frameshift_dn_ds(self):
        """Test dN and dS of codon alignments with frameshifts."""
        from Bio.codonalign.codonseq import cal_dn_ds

        codon_aln1 = codonalign.build(self.aln1, self.seqlist1)
        codon_aln2 = codonalign.build(self.aln2, self.seqlist2)
        # build pads with gaps the codons that a frameshift cuts short
        self.assertEqual(codon_aln1[1].seq[21:24], "CG-")
        self.assertEqual(codon_aln2[1].seq[15:18], "TT-")
        self.assertEqual(codon_aln2[1].seq[99:102], "AA-")
        self.assertEqual(codon_aln2[2].seq[48:51], "C--")
        dN, dS = cal_dn_ds(codon_aln1[0], codon_aln1[1], method="NG86")
        self.assertAlmostEqual(dN, 0.0463, places=4)
        self.assertAlmostEqual(dS, 0.1544, places=4)
        # The sequences differ only at the frameshifts, which are skipped
        dn, ds = codon_aln2.get_dn_ds_matrix(method="LWL85")
        self.assertEqual(dn.matrix, [[0], [0, 0], [0, 0, 0]])
        self.assertEqual(ds.matrix, [[0], [0, 0], [0, 0, 0]])

    def test_build_insertion_frameshift_dn_ds(self):
        """Test dN and dS of a frameshift insertion against frameshift deletions."""
        seq3, seq4, seq5 = self.seqlist2
        pro = SeqRecord(self.aln2[0].seq, id="pro4")
        # A copy of seq3 with one nucleotide inserted, which build pads
        # with gaps, adding gaps to the other rows against it
        inserted = SeqRecord(seq3.seq[:150] + "T" + seq3.seq[150:], id="pro4")
        codon_aln = codonalign.build(
            MultipleSeqAlignment([*self.aln2, pro]), [seq3, seq4, seq5, inserted]
        )
        self.assertEqual(codon_aln[3].seq[150:153], "T--")
        # The added gaps leave the codons of the frameshift deletions where
        # they were. This lost them, and raised ValueError.
        for rec, plain in zip(codon_aln, codonalign.build(self.aln2, self.seqlist2)):
            self.assertEqual(rec.seq[150:153], "---")
            self.assertEqual(rec.seq.rf_table, plain.seq.rf_table)
        zeros = [[0], [0, 0], [0, 0, 0], [0, 0, 0, 0]]
        for method in ("NG86", "LWL85"):
            with self.subTest(method=method):
                dn, ds = codon_aln.get_dn_ds_matrix(method=method)
                self.assertEqual(dn.matrix, zeros)
                self.assertEqual(ds.matrix, zeros)
        # Six nucleotides inserted in frame: no codon starts in their
        # columns, so they are not compared
        inserted = SeqRecord(seq3.seq[:60] + "ACGTAC" + seq3.seq[60:], id="pro4")
        codon_aln = codonalign.build(
            MultipleSeqAlignment([self.aln2[0], pro]), [seq3, inserted]
        )
        self.assertEqual(codon_aln[0].seq[60:66], "------")
        self.assertEqual(codon_aln[1].seq[60:66], "ACGTAC")
        dn, ds = codon_aln.get_dn_ds_matrix()
        self.assertEqual(dn.matrix, [[0], [0, 0]])
        self.assertEqual(ds.matrix, [[0], [0, 0]])


class Test_dn_ds(unittest.TestCase):
    def setUp(self):
        nucl = SeqIO.parse(TEST_ALIGN_FILE6[0][0], "fasta")
        prot = AlignIO.read(TEST_ALIGN_FILE6[0][1], "clustal")
        with open(TEST_ALIGN_FILE6[0][2]) as handle:
            id_corr = {i.split()[0]: i.split()[1] for i in handle}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            aln = codonalign.build(prot, nucl, corr_dict=id_corr)
        self.aln = aln

    def test_dn_ds(self):
        from Bio.codonalign.codonseq import cal_dn_ds

        codon_seq1 = self.aln[0]
        codon_seq2 = self.aln[1]
        dN, dS = cal_dn_ds(codon_seq1, codon_seq2, method="NG86")
        self.assertAlmostEqual(dN, 0.0209, places=4)
        self.assertAlmostEqual(dS, 0.0178, places=4)
        dN, dS = cal_dn_ds(codon_seq1, codon_seq2, method="LWL85")
        self.assertAlmostEqual(dN, 0.0211, places=4)
        self.assertAlmostEqual(dS, 0.0170, places=4)

        try:
            import scipy
        except ImportError:
            # Silently skip the rest of the test
            return

        # This should be present:
        from scipy.linalg import expm

        dN, dS = cal_dn_ds(codon_seq1, codon_seq2, method="YN00")
        # as PAML's yn00
        self.assertAlmostEqual(dN, 0.0201, places=4)
        self.assertAlmostEqual(dS, 0.0208, places=4)

        try:
            # New in scipy v0.11
            from scipy.optimize import minimize

            dN, dS = cal_dn_ds(codon_seq1, codon_seq2, method="ML")
            self.assertAlmostEqual(dN, 0.0194, places=4)
            self.assertAlmostEqual(dS, 0.0217, places=4)
        except ImportError:
            # TODO - Show a warning?
            pass

    def test_dn_ds_rejects_u(self):
        """Test that NG86 rejects a U rather than count it as differing from T."""
        from Bio.codonalign.codonseq import cal_dn_ds

        seq1 = codonalign.CodonSeq("AUGAAACCC")
        seq2 = codonalign.CodonSeq("ATGAAGCCC")
        with self.assertRaises(RuntimeError):
            cal_dn_ds(seq1, seq2, method="NG86")

    def test_dn_ds_matrix(self):
        # NG86 method with default codon table
        dn_correct = [
            0,
            0.02090783050583131,
            0,
            0.6115239249238438,
            0.6102203266798018,
            0,
            0.6140350835631757,
            0.6040168621204747,
            0.041180350405913294,
            0,
            0.6141532531400524,
            0.6018263135601294,
            0.06701051445629494,
            0.061470360954086874,
            0,
            0.6187088340904762,
            0.6068687248870475,
            0.07386903034833081,
            0.07357890927918581,
            0.05179847072570129,
            0,
        ]
        ds_correct = [
            0,
            0.01783718763890243,
            0,
            2.9382055377913687,
            3.0375115405379267,
            0,
            2.008913071877126,
            2.0182088023715616,
            0.5638033197005285,
            0,
            2.771425931736778,
            2.7353083173058295,
            0.6374483799734671,
            0.723542095485497,
            0,
            # Too many synonymous differences between the rat and the
            # fly sequences for NG86 to correct; this was -1
            math.nan,
            math.nan,
            0.953865978141643,
            1.182154857347706,
            0.843182957978177,
            0,
        ]
        with self.assertWarnsRegex(BiopythonWarning, "^NG86 cannot"):
            dn, ds = self.aln.get_dn_ds_matrix()
        dn_list = []
        for i in dn.matrix:
            dn_list.extend(i)
        for dn_cal, dn_corr in zip(dn_list, dn_correct):
            self.assertAlmostEqual(dn_cal, dn_corr, places=4)
        ds_list = []
        for i in ds.matrix:
            ds_list.extend(i)
        for ds_cal, ds_corr in zip(ds_list, ds_correct):
            if math.isnan(ds_corr):
                self.assertTrue(math.isnan(ds_cal))
            else:
                self.assertAlmostEqual(ds_cal, ds_corr, places=4)
        # LWL85 method with user specified codon table
        dn_correct = [
            0,
            0.02040017274543761,
            0,
            0.6357565704110552,
            0.6350555099061804,
            0,
            0.6326332180170293,
            0.6235817405000689,
            0.046368959735469464,
            0,
            0.6475034783787625,
            0.6372033288345671,
            0.07344045228462681,
            0.07223785142164467,
            0,
            0.6399664487237698,
            0.6278036887068891,
            0.08310311966390273,
            0.08568698005598897,
            0.059944126284334445,
            0,
        ]
        ds_correct = [
            0,
            0.019580765877653315,
            0,
            2.031857057230179,
            2.0668216131397137,
            0,
            1.699660630564657,
            1.698605472271965,
            0.5209710723015657,
            0,
            1.8622121922473729,
            1.8344006980089567,
            0.5748202364140963,
            0.6201527130491078,
            0,
            2.3594491316933413,
            2.4291632003286447,
            0.8077681072185632,
            0.9491264436707163,
            0.7334848415191596,
            0,
        ]
        dn, ds = self.aln.get_dn_ds_matrix(
            method="LWL85", codon_table=CodonTable.unambiguous_dna_by_id[3]
        )
        dn_list = []
        for i in dn.matrix:
            dn_list.extend(i)
        for dn_cal, dn_corr in zip(dn_list, dn_correct):
            self.assertAlmostEqual(dn_cal, dn_corr, places=4)
        ds_list = []
        for i in ds.matrix:
            ds_list.extend(i)
        for ds_cal, ds_corr in zip(ds_list, ds_correct):
            self.assertAlmostEqual(ds_cal, ds_corr, places=4)

    def test_dn_ds_matrix_yn00(self):
        """Test get_dn_ds_matrix with the YN00 method."""
        try:
            import scipy
        except ImportError:
            self.skipTest("the YN00 method needs scipy")
        # YN00 method: the values are those of the yn00 program of PAML
        # 4.10.10, run with weighting = 1 on each pair. yn00 prints four decimal
        # places and stops iterating at a coarser tolerance, so they agree to
        # within 1e-4. This raised "ValueError: math domain error".
        dn_yn00 = [
            [0.0201],
            [0.5653, 0.5633],
            [0.5554, 0.5443, 0.0389],
            [0.5616, 0.5480, 0.0653, 0.0589],
            [0.5733, 0.5591, 0.0718, 0.0712, 0.0513],
        ]
        ds_yn00 = [
            [0.0208],
            [4.9068, 4.9042],
            [4.8692, 4.8612, 0.6970],
            [4.9004, 4.8908, 0.6427, 0.8678],
            [4.9226, 4.9138, 0.9766, 2.6007, 0.8877],
        ]
        dn, ds = self.aln.get_dn_ds_matrix(method="YN00")
        for i, (dn_row, ds_row) in enumerate(zip(dn_yn00, ds_yn00), 1):
            for j, (dn_corr, ds_corr) in enumerate(zip(dn_row, ds_row)):
                self.assertAlmostEqual(dn[i, j], dn_corr, delta=1e-4)
                self.assertAlmostEqual(ds[i, j], ds_corr, delta=1e-4)

    def test_dn_ds_tree(self):
        """Test get_dn_ds_tree builds UPGMA trees by default, or WPGMA trees."""
        from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

        constructor = DistanceTreeConstructor()
        dn, ds = self.aln.get_dn_ds_matrix(method="LWL85")
        dn_tree, ds_tree = self.aln.get_dn_ds_tree(dn_ds_method="LWL85")
        # UPGMA joins the fly sequences 0 and 1 to the vertebrate sequences
        # 2 to 5 at the mean of the eight dS distances between them, counting
        # every pair equally.
        mean = sum(ds[i, j] for i in (0, 1) for j in (2, 3, 4, 5)) / 8
        for i in (0, 1):
            for j in (2, 3, 4, 5):
                with self.subTest(i=i, j=j):
                    distance = ds_tree.distance(ds.names[i], ds.names[j])
                    self.assertAlmostEqual(distance, mean)
        for tree, dm in ((dn_tree, dn), (ds_tree, ds)):
            self.assertEqual(
                format(tree, "newick"), format(constructor.upgma(dm), "newick")
            )
        dn_tree, ds_tree = self.aln.get_dn_ds_tree(
            dn_ds_method="LWL85", tree_method="WPGMA"
        )
        for tree, dm in ((dn_tree, dn), (ds_tree, ds)):
            self.assertEqual(
                format(tree, "newick"), format(constructor.wpgma(dm), "newick")
            )
        # NG86 cannot give dS between the rat and the fly sequences. With -1
        # for it, UPGMA joined the rat to the flies.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            with self.assertRaisesRegex(
                ValueError, r"RATEGFR and gi\|24657088\|ref\|NM_057410\.3\| \(dS\)"
            ):
                self.aln.get_dn_ds_tree()


class Test_dn_ds_too_divergent(unittest.TestCase):
    """Test dN and dS of pairs too divergent for NG86 or LWL85."""

    def test_cal_dn_ds(self):
        """Test the two pairs of biopython/biopython#1611."""
        from Bio.codonalign.codonseq import cal_dn_ds

        # NG86 gave -1 for both, and LWL85 raised a math domain error
        seq1 = codonalign.CodonSeq("TTTTTTTTT")
        seq2 = codonalign.CodonSeq("CCACCACCA")
        for method in ("NG86", "LWL85"):
            with self.subTest(method=method):
                with self.assertWarnsRegex(BiopythonWarning, f"^{method} cannot"):
                    dn, ds = cal_dn_ds(seq1, seq2, method=method)
                self.assertTrue(math.isnan(dn))
                self.assertTrue(math.isnan(ds))
        # Only the twofold degenerate sites have too many differences
        seq1 = codonalign.CodonSeq("ATGGGTACCGTGTTGTAT")
        seq2 = codonalign.CodonSeq("ATGGCTACCGTGTTAGTA")
        with self.assertWarnsRegex(BiopythonWarning, "returning nan for dS$"):
            dn, ds = cal_dn_ds(seq1, seq2, method="LWL85")
        self.assertAlmostEqual(dn, 0.4257, places=4)
        self.assertTrue(math.isnan(ds))

    def setUp(self):
        self.aln = codonalign.CodonAlignment(
            [
                SeqRecord(codonalign.CodonSeq("TTTGGGGGGGGG"), id="a"),
                SeqRecord(codonalign.CodonSeq("TTCGGGGGGGGG"), id="b"),
                SeqRecord(codonalign.CodonSeq("TTTTTTTTTTTT"), id="c"),
            ]
        )

    def test_dn_ds_matrix(self):
        """Test the matrices hold nan for the undefined distances only."""
        with self.assertWarns(BiopythonWarning):
            dn, ds = self.aln.get_dn_ds_matrix()
        self.assertEqual(dn["a", "b"], 0)
        self.assertGreater(ds["a", "b"], 0)
        self.assertTrue(math.isnan(dn["a", "c"]))
        self.assertGreater(ds["a", "c"], 0)
        self.assertTrue(math.isnan(dn["b", "c"]))
        self.assertTrue(math.isnan(ds["b", "c"]))

    def test_dn_ds_tree(self):
        """Test get_dn_ds_tree refuses undefined distances, naming them.

        With -1 for them, UPGMA joined c to a, the pair furthest apart.
        """
        message = (
            r"^NG86 gives no distance for c and a \(dN\), c and b \(dN\), "
            r"c and b \(dS\), as they are too divergent"
        )
        for tree_method in ("UPGMA", "WPGMA", "NJ"):
            with self.subTest(tree_method=tree_method):
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", BiopythonWarning)
                    with self.assertRaisesRegex(ValueError, message):
                        self.aln.get_dn_ds_tree(tree_method=tree_method)


class Test_dn_ds_not_codon_aligned(unittest.TestCase):
    """Test dN and dS of alignments whose gaps split codons."""

    def setUp(self):
        # Test set 3 of Test_build
        self.codon_table = CodonTable.unambiguous_dna_by_id[3]
        self.seq1 = "ATGGCAAGGGACCACCCAGTTGGGCACTGATATGATCGGGTGTATTTGCAGAGTAGTAACCTTTCTTTTCTCAAGACCATCCAG"
        self.seq2 = "ATGGCAAGGCACCATCCAGTTGAGCACTGATATGATCGGGTGTATTTGCAGAGTAGTAACGTGTCTCTGCTCAAGACCATCCAG"

    def assertNotCodonAligned(self, seq1, seq2, columns):
        """Assert that cal_dn_ds rejects the pair, naming the columns."""
        from Bio.codonalign.codonseq import cal_dn_ds

        seq1 = codonalign.CodonSeq(seq1)
        seq2 = codonalign.CodonSeq(seq2)
        message = f"^Not a codon alignment: alignment columns {columns} "
        for method in ("NG86", "LWL85", "YN00", "ML"):
            with self.subTest(method=method):
                with self.assertRaisesRegex(ValueError, message):
                    cal_dn_ds(seq1, seq2, method=method)
                with self.assertRaisesRegex(ValueError, message):
                    cal_dn_ds(seq2, seq1, method=method)

    def test_out_of_frame(self):
        """Test codons put out of frame by a gap that splits a codon."""
        from Bio.codonalign.codonseq import cal_dn_ds

        # The codons of the first are AAA, CCC, GGG and TTT, but the
        # three columns 10-12 hold GTT. This gave dN of 0.165 by LWL85 and
        # 0.135 by ML for two sequences with no differences, and NG86 and
        # YN00 raised TypeError and ValueError.
        self.assertNotCodonAligned("AAACCC-GGGTTT--", "AAACCCGGGTTT---", "10-12")
        # This raised IndexError
        self.assertNotCodonAligned("AAACCC--GGGTT----T", "AAACCCGGGTTT------", "10-12")
        # The arguments are checked first
        with self.assertRaisesRegex(RuntimeError, "^cfreq can only be specified"):
            cal_dn_ds(
                codonalign.CodonSeq("AAACCC-GGGTTT--"),
                codonalign.CodonSeq("AAACCCGGGTTT---"),
                cfreq="F61",
            )

    def test_issue_pair(self):
        """Test the nucleotide alignment of biopython/biopython#4999."""
        # The 62 gaps after the first 140 nucleotides split a codon, and put
        # the codons after them out of frame. This raised "ValueError:
        # Sequence length is not a multiple of three", though it is.
        seq1 = (
            "-" * 9
            + (
                "atgaccgccgaatcgccgacgatccgcctcgaacgctacagcgaacgccacgtggaaggcctca"
                "ccgccctgtacaacgaccccgccgtggcccgccaggtcctgcagatgccctaccagtcggtgga"
                "gcagcgccgcaa"
            )
            + "-" * 62
            + (
                "gcgcctgcacgacagcgcggacgacgaccgcctgctgatcctcgtcgccctgcaccagggggac"
                "gtcatcggcagcgccagcctcgaacagcacccgcgcatccgccg"
            )
            + "-" * 6
            + (
                "cagccacagcggcagcatcggcatgggcgtcgccgtcgcctggcagggcaagggcgtcggcagc"
                "cgcctgctcggcgaactgctggacatcgccgacaactggatgaacctgcgccgcgtcgaactca"
                "ccgtctacaccgacaacgccccggcgctggcgctgtatcgcaagttcggcttcgagaccgaagg"
                "agagatgcgcgactacgcggtgcgggatgggcgcttcgtggacgtgtacagcatggcgcgcttg"
                "cggag"
            )
            + "-" * 16
            + "ggtcgaggggcgg"
            + "-" * 79
            + "gtgg"
            + "-" * 1
            + "gggagtga"
            + "-" * 181
        )
        # The other sequence differs from it in three columns
        seq2 = seq1[:68] + "t" + seq1[69:347] + "c" + seq1[348:391] + "t" + seq1[392:]
        self.assertEqual(len(seq1), 888)
        self.assertNotCodonAligned(seq1, seq2, "211-213")

    def test_split_codon_skipped(self):
        """Test a codon split by gaps that leave the codons after it in frame.

        It is skipped, like any codon with a gap.
        """
        from Bio.codonalign.codonseq import cal_dn_ds

        seq1, seq2, codon_table = self.seq1, self.seq2, self.codon_table
        # Split the codon at i as XY- Z--, with the other sequence's codon
        # against XY-, or remove it, leaving six gaps. With the codon at 81,
        # the last one, the split raised IndexError.
        for i in (30, 81):
            x, y, z = seq1[i : i + 3]
            split = codonalign.CodonSeq(f"{seq1[:i]}{x}{y}-{z}--{seq1[i + 3 :]}")
            gapped = codonalign.CodonSeq(f"{seq1[:i]}------{seq1[i + 3 :]}")
            other = codonalign.CodonSeq(f"{seq2[: i + 3]}---{seq2[i + 3 :]}")
            for method in ("NG86", "LWL85"):
                with self.subTest(i=i, method=method):
                    self.assertEqual(
                        cal_dn_ds(split, other, method, codon_table),
                        cal_dn_ds(gapped, other, method, codon_table),
                    )

    def test_no_codon_start_skipped(self):
        """Test three columns in which no codon of the reading frame starts.

        They are skipped, like a codon with a gap.
        """
        from Bio.codonalign.codonseq import cal_dn_ds

        seq1, codon_table = self.seq1, self.codon_table
        other = codonalign.CodonSeq(self.seq2)
        # A reading frame that skips the codon at i, as for an insertion.
        # This raised TypeError, or IndexError for the last codon, at 81.
        for i in (30, 81):
            rf_table = [j for j in range(0, len(seq1), 3) if j != i]
            skipped = codonalign.CodonSeq(seq1, rf_table=rf_table)
            gapped = codonalign.CodonSeq(f"{seq1[:i]}---{seq1[i + 3 :]}")
            for method in ("NG86", "LWL85"):
                with self.subTest(i=i, method=method):
                    self.assertEqual(
                        cal_dn_ds(skipped, other, method, codon_table),
                        cal_dn_ds(gapped, other, method, codon_table),
                    )

    def test_alignment_methods(self):
        """Test the CodonAlignment methods and mktest on codons out of frame."""
        codon_aln = codonalign.CodonAlignment(
            [
                SeqRecord(codonalign.CodonSeq("AAACCC-GGGTTT--"), id="alpha"),
                SeqRecord(codonalign.CodonSeq("AAACCCGGGTTT---"), id="beta"),
                SeqRecord(codonalign.CodonSeq("AAACCAGGGTTT---"), id="gamma"),
            ]
        )
        # These raised TypeError, and mktest gave a p-value
        message = "^Not a codon alignment: alignment columns 10-12 "
        with self.assertRaisesRegex(ValueError, message) as cm:
            codon_aln.get_dn_ds_matrix()
        self.assertIn(" Bio.Align.CodonAligner makes", str(cm.exception))
        self.assertIn(" Bio.Align.analysis.calculate_dn_ds ", str(cm.exception))
        with self.assertRaisesRegex(ValueError, message):
            codon_aln.get_dn_ds_tree()
        with self.assertRaisesRegex(ValueError, message):
            codonalign.mktest([codon_aln, codon_aln])


try:
    import numpy as np
except ImportError:
    np = None

if np:

    class Test_MK(unittest.TestCase):
        def test_mk(self):
            p = SeqIO.index(TEST_ALIGN_FILE7[0][0], "fasta")
            pro_aln = AlignIO.read(TEST_ALIGN_FILE7[0][1], "clustal")
            codon_aln = codonalign.build(pro_aln, p)
            p.close()  # Close indexed FASTA file
            # Same Drosophila adh data, and the same value, as
            # test_Align_codonalign.py asserts for Bio.Align.analysis.mktest.
            # The two implementations agreeing to every digit is the point:
            # they were fixed separately from the same reading of the data.
            # The value before the _prim fix was 0.0020645725725430097, from a
            # contingency table inflated by an extra substitution at every
            # site with three or more distinct codons.
            self.assertAlmostEqual(
                codonalign.mktest([codon_aln[1:12], codon_aln[12:16], codon_aln[16:]]),
                0.0023196627124160174,
            )

        def test_mk_skips_gapped_codons(self):
            """Test that a codon position with a gap in any sequence is skipped."""

            def mktest(*species):
                return codonalign.mktest(
                    [
                        codonalign.CodonAlignment(
                            SeqRecord(codonalign.CodonSeq(row), id=f"{i}_{j}")
                            for j, row in enumerate(rows)
                        )
                        for i, rows in enumerate(species)
                    ]
                )

            # Six codon positions: two synonymous polymorphisms, then one
            # synonymous and three nonsynonymous fixed differences.
            species1 = ["CTTGCTAAAGATTTTATG", "CTCGCCAAAGATTTTATG"]
            species2 = ["CTTGCTAAGGAATTAACG", "CTTGCTAAGGAATTAACG"]
            pvalue = mktest(species1, species2)
            # G test of the table [1, 3, 2, 0], worked out by hand
            g = math.log(1 / 2) + 3 * math.log(3 / 2) + 2 * math.log(2)
            self.assertAlmostEqual(pvalue, math.erfc(math.sqrt(g)))
            # One more position, TGG against TGT, would add a fixed difference.
            # With a gap in one sequence it is left out instead; it used to
            # raise KeyError.
            self.assertNotEqual(
                mktest(
                    ["TGG" + row for row in species1], ["TGT" + row for row in species2]
                ),
                pvalue,
            )
            self.assertEqual(
                mktest(
                    ["---" + species1[0], "TGG" + species1[1]],
                    ["TGT" + row for row in species2],
                ),
                pvalue,
            )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
