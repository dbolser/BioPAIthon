# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Malformed input to the Bio.Align parsers and writers raises ValueError.

These checks used to be ``assert`` statements, which vanish under
``python -O``, so malformed input then gave silently wrong results; with
asserts enabled, the bare ``AssertionError`` did not say what was wrong.
Each test feeds one malformed input and expects a ValueError naming it.
"""

import tempfile
import unittest
from io import BytesIO
from io import StringIO

import numpy as np

import support

from Bio import Align
from Bio.Align import Alignment
from Bio.Align import Alignments
from Bio.Align import bigbed
from Bio.Align import substitution_matrices
from Bio.Align.analysis import _count_site_YN00
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse(text, fmt):
    """Parse all alignments in text, given as str or bytes."""
    if isinstance(text, bytes):
        stream = BytesIO(text)
    else:
        stream = StringIO(text)
    return list(Align.parse(stream, fmt))


def replace_once(text, old, new):
    """Replace old by new in text, which must contain old exactly once."""
    if text.count(old) != 1:
        raise AssertionError(f"expected exactly one {old!r} in the test input")
    return text.replace(old, new)


class TestA2M(unittest.TestCase):
    def test_upper_case_letter_in_insert_column(self):
        """A match letter where the first sequence has an insertion."""
        text = ">seq1\nACDef\n>seq2\nACDEf\n"
        with self.assertRaises(ValueError) as cm:
            parse(text, "a2m")
        self.assertIn("column 4 of sequence seq2, found 'E'", str(cm.exception))


class TestAnalysis(unittest.TestCase):
    def test_count_site_yn00_unequal_codon_counts(self):
        """YN00 site counting needs as many codons in both sequences."""
        codon_table = CodonTable.generic_by_id[1]
        with self.assertRaises(ValueError) as cm:
            _count_site_YN00(["ATG", "AAA"], ["ATG"], {}, 1, codon_table)
        self.assertIn("found 2 and 1", str(cm.exception))


class TestBigBed(unittest.TestCase):
    def test_autosql_field_without_semicolon(self):
        """Every AutoSQL field definition ends with a semicolon."""
        text = (
            'table bed\n"Browser Extensible Data"\n(\n string chrom "Chromosome"\n)\n'
        )
        with self.assertRaises(ValueError) as cm:
            bigbed.AutoSQLTable.from_string(text)
        self.assertIn("'string chrom' to end with ';'", str(cm.exception))

    def test_autosql_unknown_data_type(self):
        text = 'table bed\n"Browser Extensible Data"\n(\n text chrom; "Chromosome"\n)\n'
        with self.assertRaises(ValueError) as cm:
            bigbed.AutoSQLTable.from_string(text)
        self.assertIn("Unknown AutoSQL data type 'text'", str(cm.exception))

    def test_unsupported_version(self):
        """The header gives the bigBed format version."""
        with open(support.DATA / "Blat" / "bigbedtest.bb", "rb") as stream:
            data = stream.read()
        # bytes 4-5 hold the version number; bigbedtest.bb is little-endian
        data = data[:4] + (3).to_bytes(2, "little") + data[6:]
        with self.assertRaises(ValueError) as cm:
            parse(data, "bigbed")
        self.assertIn("Expected bigBed version 4, found 3", str(cm.exception))

    def test_item_without_nul_terminator(self):
        """A data item that does not end in a NUL byte.

        With one item per slot, searching hands each whole data block to the
        item parser, which would silently drop its last byte if the block was
        not NUL-terminated.
        """
        target = SeqRecord(Seq(None, length=1000), id="chr1")
        query = SeqRecord(Seq(None, length=50), id="uniquename")
        alignment = Alignment([target, query], np.array([[100, 150], [0, 50]]))
        alignments = Alignments([alignment])
        alignments.targets = [target]
        with tempfile.TemporaryFile() as stream:
            Align.write(
                alignments, stream, "bigbed", bedN=4, compress=False, itemsPerSlot=1
            )
            stream.seek(0)
            data = stream.read()
        data = replace_once(data, b"uniquename\0", b"uniquename!")
        alignments = Align.parse(BytesIO(data), "bigbed")
        with self.assertRaises(ValueError) as cm:
            list(alignments.search("chr1", 100, 150))
        self.assertIn("to end with a NUL byte, found b'!'", str(cm.exception))


class TestBigMaf(unittest.TestCase):
    def test_writing_reference_on_reverse_strand(self):
        """The reference must be aligned on its forward strand."""
        reference = SeqRecord(Seq("ACGTACGTAC"), id="hg38.chr1")
        other = SeqRecord(Seq("ACGTACGTAC"), id="mm10.chr2")
        alignment = Alignment([reference, other], np.array([[10, 0], [0, 10]]))
        alignments = Alignments([alignment])
        alignments.targets = [reference]
        with self.assertRaises(ValueError) as cm:
            Align.write(alignments, BytesIO(), "bigmaf")
        self.assertIn("hg38.chr1 to be aligned on the forward", str(cm.exception))


class TestBigPsl(unittest.TestCase):
    def test_writing_mixed_step_sizes(self):
        """Blocks are all nucleotide-nucleotide or all translated."""
        target = SeqRecord(Seq(None, length=100), id="chr1")
        query = SeqRecord(Seq(None, length=10), id="query")
        coordinates = np.array([[0, 6, 10], [0, 2, 6]])
        alignment = Alignment([target, query], coordinates)
        alignments = Alignments([alignment])
        alignments.targets = [target]
        with self.assertRaises(ValueError) as cm:
            Align.write(alignments, BytesIO(), "bigpsl")
        self.assertIn("found steps 4, 4", str(cm.exception))


class TestClustal(unittest.TestCase):
    def test_identifier_changes_between_blocks(self):
        text = """\
CLUSTAL W (1.83) multiple sequence alignment


seq1      ACGT
seq2      ACGT

seq1      ACGT
seq3      ACGT
"""
        with self.assertRaises(ValueError) as cm:
            parse(text, "clustal")
        self.assertIn("Expected sequence seq2", str(cm.exception))

    def test_unequal_lengths_in_first_block(self):
        text = """\
CLUSTAL W (1.83) multiple sequence alignment


seq1      ACGT
seq2      ACG

"""
        with self.assertRaises(ValueError) as cm:
            parse(text, "clustal")
        self.assertIn("Expected 4 columns for seq2", str(cm.exception))


class TestEmboss(unittest.TestCase):
    def test_end_coordinate_does_not_match_letters(self):
        with open(support.DATA / "Emboss" / "needle.txt") as stream:
            text = stream.read()
        text = replace_once(
            text,
            "ref_rec            1 KILIVDD----QYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDL     46",
            "ref_rec            1 KILIVDD----QYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDL     47",
        )
        with self.assertRaises(ValueError) as cm:
            parse(text, "emboss")
        self.assertIn(
            "coordinates of ref_rec to span the 46 letters", str(cm.exception)
        )


class TestExonerate(unittest.TestCase):
    header = "Command line: [exonerate]\nHostname: [host]\n"

    def test_missing_command_line(self):
        text = "Hostname: [host]\n-- completed exonerate analysis\n"
        with self.assertRaises(ValueError) as cm:
            parse(text, "exonerate")
        self.assertIn("to start with 'Command line: '", str(cm.exception))

    def test_vulgar_splice_site_length(self):
        """Splice sites are two nucleotides long."""
        text = (
            self.header
            + "vulgar: q 0 10 + t 0 15 + 30 M 5 5 5 0 3 I 0 2 3 0 2 M 5 5\n"
            + "-- completed exonerate analysis\n"
        )
        with self.assertRaises(ValueError) as cm:
            parse(text, "exonerate")
        self.assertIn("Expected a length of 2 for operation '5'", str(cm.exception))

    def test_writing_deletion_with_query_step(self):
        """A deletion consumes the target only."""
        target = SeqRecord(Seq(None, length=20), id="t")
        query = SeqRecord(Seq(None, length=20), id="q")
        alignment = Alignment([target, query], np.array([[0, 5, 10], [0, 5, 10]]))
        alignment.score = 10
        alignment.operations = bytearray(b"MD")
        with self.assertRaises(ValueError) as cm:
            format(alignment, "exonerate")
        self.assertIn("query step of 0 for operation 'D'", str(cm.exception))


class TestMaf(unittest.TestCase):
    header = "##maf version=1\n\n"

    def test_i_line_for_another_sequence(self):
        text = (
            self.header
            + "a score=0\n"
            + "s hg38.chr1 10 4 + 100 ACGT\n"
            + "s mm10.chr2 20 4 + 100 ACGT\n"
            + "i hg38.chr1 C 0 C 0\n"
        )
        with self.assertRaises(ValueError) as cm:
            parse(text, "maf")
        self.assertIn("Expected 'i' line for mm10.chr2", str(cm.exception))
        self.assertIn("i hg38.chr1 C 0 C 0", str(cm.exception))

    def test_e_line_with_unknown_status(self):
        text = (
            self.header
            + "a score=0\n"
            + "s hg38.chr1 10 4 + 100 ACGT\n"
            + "e mm10.chr2 20 4 + 100 X\n"
        )
        with self.assertRaises(ValueError) as cm:
            parse(text, "maf")
        self.assertIn("in 'e' line, found 'X'", str(cm.exception))

    def test_first_block_without_a_line(self):
        text = self.header + "s hg38.chr1 10 4 + 100 ACGT\n"
        with self.assertRaises(ValueError) as cm:
            parse(text, "maf")
        self.assertIn("Expected an 'a' line", str(cm.exception))


class TestMauve(unittest.TestCase):
    def test_sequence_numbers_out_of_order(self):
        text = """\
#FormatVersion Mauve1
#Sequence1File	a.fa
#Sequence3File	b.fa
> 1:1-4 + a.fa
ACGT
=
"""
        with self.assertRaises(ValueError) as cm:
            parse(text, "mauve")
        self.assertIn("Expected #Sequence2File", str(cm.exception))

    def test_unknown_strand(self):
        text = """\
#FormatVersion Mauve1
#Sequence1File	a.fa
#Sequence2File	b.fa
> 1:1-4 * a.fa
ACGT
=
"""
        with self.assertRaises(ValueError) as cm:
            parse(text, "mauve")
        self.assertIn("Expected strand '+' or '-', found '*'", str(cm.exception))


class TestMsf(unittest.TestCase):
    def test_text_after_sequence_data(self):
        with open(support.DATA / "msf" / "DOA_prot.msf") as stream:
            text = stream.read()
        text += "\nsome trailing text\n"
        with self.assertRaises(ValueError) as cm:
            parse(text, "msf")
        self.assertIn("found 'some trailing text\\n'", str(cm.exception))


class TestPhylip(unittest.TestCase):
    def test_interleaved_block_missing_a_line(self):
        text = """\
3 12
seq1      ACGT
seq2      ACGT
seq3      ACGT

ACGT
ACGT

ACGT
ACGT
ACGT
"""
        with self.assertRaises(ValueError) as cm:
            parse(text, "phylip")
        self.assertIn(
            "Expected 3 sequence lines in each block of the interleaved "
            "alignment, found 2",
            str(cm.exception),
        )


class TestPsl(unittest.TestCase):
    def test_writing_mixed_step_sizes(self):
        """Blocks are all nucleotide-nucleotide or all translated."""
        target = SeqRecord(Seq(None, length=100), id="chr1")
        query = SeqRecord(Seq(None, length=10), id="query")
        coordinates = np.array([[0, 2, 3], [0, 1, 3]])
        alignment = Alignment([target, query], coordinates)
        with self.assertRaises(ValueError) as cm:
            format(alignment, "psl")
        self.assertIn("found steps 2, 1", str(cm.exception))


class TestSam(unittest.TestCase):
    header = "@SQ\tSN:chr1\tLN:100\n"

    def test_score_tag_of_wrong_type(self):
        """The AS tag holds an integer."""
        text = self.header + "read1\t0\tchr1\t1\t255\t4M\t*\t0\t0\tACGT\t*\tAS:Z:4\n"
        with self.assertRaises(ValueError) as cm:
            parse(text, "sam")
        self.assertIn("Expected type 'i' for tag AS of read1", str(cm.exception))

    def test_sequence_shorter_than_cigar(self):
        text = self.header + "read1\t0\tchr1\t1\t255\t5M\t*\t0\t0\tACGT\t*\n"
        with self.assertRaises(ValueError) as cm:
            parse(text, "sam")
        self.assertIn("Expected 5 letters in the sequence of read1", str(cm.exception))

    def test_header_tag_not_two_letters(self):
        text = "@HD\tVN:1.6\tSORT:coordinate\n"
        with self.assertRaises(ValueError) as cm:
            parse(text, "sam")
        self.assertIn("found 'SORT'", str(cm.exception))


class TestStockholm(unittest.TestCase):
    def test_gr_line_for_another_sequence(self):
        text = """\
# STOCKHOLM 1.0
seq1 ACGT
seq2 ACGT
#=GR seq1 SS ....
//
"""
        with self.assertRaises(ValueError) as cm:
            parse(text, "stockholm")
        self.assertIn("Expected #=GR line for seq2", str(cm.exception))

    def test_gap_characters_disagree_in_a_column(self):
        """'.' marks insert columns and '-' match columns, never both."""
        text = """\
# STOCKHOLM 1.0
seq1 AC.T
seq2 AC-T
//
"""
        with self.assertRaises(ValueError) as cm:
            parse(text, "stockholm")
        self.assertIn("column 3 of seq2", str(cm.exception))

    def test_writing_operations_of_wrong_length(self):
        alignment = Alignment(
            [SeqRecord(Seq("ACGT"), id="seq1"), SeqRecord(Seq("ACGT"), id="seq2")]
        )
        alignment.operations = bytearray(b"MMM")
        with self.assertRaises(ValueError) as cm:
            format(alignment, "stockholm")
        self.assertIn("(4), found 3", str(cm.exception))


class TestSubstitutionMatrices(unittest.TestCase):
    def test_rows_out_of_order(self):
        text = "   A  C\nC  1  0\nA  0  1\n"
        with self.assertRaises(ValueError) as cm:
            substitution_matrices.read(StringIO(text))
        self.assertIn(
            "Expected the row for 'A', in the order of the column header",
            str(cm.exception),
        )


class TestTabular(unittest.TestCase):
    def test_row_with_missing_column(self):
        text = """\
# BLASTN 2.13.0+
# Query: query1
# Database: db
# Fields: query id, subject id, % identity
# 1 hits found
query1\tsubject1
# BLAST processed 1 queries
"""
        with self.assertRaises(ValueError) as cm:
            parse(text, "tabular")
        self.assertIn("Expected 3 columns", str(cm.exception))

    def test_query_id_differs_from_header(self):
        text = """\
# BLASTN 2.13.0+
# Query: query1
# Database: db
# Fields: query id, subject id, % identity
# 1 hits found
query2\tsubject1\t100.0
# BLAST processed 1 queries
"""
        with self.assertRaises(ValueError) as cm:
            parse(text, "tabular")
        self.assertIn("Expected query id query1", str(cm.exception))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
