"""Tests for SeqIO GFA module."""

import unittest
import warnings
from io import StringIO

from Bio import BiopythonWarning
from Bio import SeqIO


class TestRead(unittest.TestCase):
    def test_read_GFA1(self):
        """Test parsing valid GFA 1.x files."""
        records = list(SeqIO.parse("GFA/seq.gfa", "gfa1"))
        self.assertEqual(len(records), 8)
        self.assertEqual(records[6].id, "MTh13014")
        self.assertEqual(
            records[6].seq,
            "TTAGGTCTCCACCCCTGACTCCCCTCAGCCATAGAAGGCCCCACCCCAGTCTCAGCCCTACTCCACTCAAGCACTATAGTTGTAGCAGGAATCTTCTTACTCATCCGCTTCCACCCCCTAGCAGAAAATAGCCCACTAATCCAAACTCTAACACTATGCTTAGGCGCTATCACCACTCTGTTCGCAGCAGTCTGCGCCCTTACACAAAATGACATCAAAAAAATCGTAGCCTTCTCCACTTCAAGTCAACTAGGACTCATAATAGTTACAATCGGCATCAACCAACCACACCTAGCATTCCTGCACATCTGTACCCACGCCTTCTTCAAAGCCATACTATTTATGTGCTCCGGGTCCATCATCCACAACCTTAACAATGAACAAGATATTCGAAAAATAGGAGGACTACTCAAAACCATACCTCTCACTTCAACCTCCCTCACCATTGGCAGCCTAGCATTAGCAGGAATACCTTTCCTCACAGGTTTCTACTCCAAAGACC",
        )
        self.assertEqual(records[0].annotations["SN"], ("Z", "MT_human"))
        self.assertEqual(records[0].annotations["SO"], ("i", "0"))

        records = list(SeqIO.parse("GFA/seq_with_len.gfa", "gfa1"))
        self.assertEqual(len(records), 9)
        self.assertEqual(
            records[8].seq,
            "GAAAAATTGCCCTTGGTTTTCGCTTCGCTCAAACTCTATTGAACTTCGCTTTCGCTCAGTTCGTCGGGGCAATTTTTTGGTTAATACTT",
        )

        records = list(SeqIO.parse("GFA/fake_with_checksum.gfa", "gfa1"))
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].seq, "AAA")

        records = list(SeqIO.parse("GFA/no_seq.gfa", "gfa1"))
        self.assertEqual(len(records), 9)
        self.assertEqual(len(records[0]), 528)

    def test_read_GFA2(self):
        """Test parsing valid GFA 2.0 files."""
        records = list(SeqIO.parse("GFA/fake_gfa2.gfa", "gfa2"))
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].seq, "AAA")

    def test_read_GFA2_no_seq(self):
        """Test a GFA 2.0 segment without sequence keeps its length."""
        record = SeqIO.read(StringIO("S\ts1\t100\t*\n"), "gfa2")
        self.assertEqual(len(record), 100)
        self.assertFalse(record.seq.defined)

    def test_read_GFA2_no_seq_LN(self):
        """Test LN gives the length of a GFA 2.0 segment without sequence."""
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            record = SeqIO.read(StringIO("S\ts1\t100\t*\tLN:i:528\n"), "gfa2")
        self.assertEqual(len(record), 528)
        self.assertFalse(record.seq.defined)

    def test_read_GFA2_no_seq_checksum(self):
        """Test a GFA 2.0 segment without sequence keeps an unverifiable checksum."""
        checksum = "ABCD" * 16
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            record = SeqIO.read(StringIO(f"S\ts1\t100\t*\tSH:H:{checksum}\n"), "gfa2")
        self.assertEqual(len(record), 100)
        self.assertEqual(record.annotations["SH"], ("H", checksum))


class TestCorrupt(unittest.TestCase):
    def test_corrupt_gfa2(self):
        """Check a GFA 1.x file does not parse in GFA 2."""
        with self.assertRaises(ValueError):
            list(SeqIO.parse("GFA/seq.gfa", "gfa2"))

    def test_corrupt_segment_fields(self):
        """Check a GFA file with invalid fields on a segment line."""
        with self.assertRaises(ValueError):
            list(SeqIO.parse("GFA/corrupt_segment_fields.gfa", "gfa1"))

    def test_corrupt_len(self):
        """Check a GFA file with an incorrect length."""
        with self.assertWarns(BiopythonWarning):
            list(SeqIO.parse("GFA/corrupt_len.gfa", "gfa1"))

    def test_corrupt_repeated_len(self):
        """Check only the first LN tag sets the length of a segment without sequence."""
        for fmt, line in (
            ("gfa1", "S\ts1\t*\tLN:i:5\tLN:i:6\n"),
            ("gfa2", "S\ts1\t100\t*\tLN:i:5\tLN:i:6\n"),
        ):
            with self.subTest(fmt=fmt):
                with self.assertWarnsRegex(BiopythonWarning, "incorrect length"):
                    record = SeqIO.read(StringIO(line), fmt)
                self.assertEqual(len(record), 5)

    def test_corrupt_negative_len_gfa2(self):
        """Check a GFA 2.0 segment without sequence cannot have a negative length."""
        with self.assertRaisesRegex(ValueError, "non-negative length: S\ts1\t-5\t\\*"):
            SeqIO.read(StringIO("S\ts1\t-5\t*\n"), "gfa2")

    def test_corrupt_checksum(self):
        """Check a GFA file with an incorrect checksum."""
        with self.assertWarns(BiopythonWarning):
            list(SeqIO.parse("GFA/corrupt_checksum.gfa", "gfa1"))

    def test_corrupt_tag_name(self):
        """Check a GFA file with an invalid tag name."""
        with self.assertWarns(BiopythonWarning):
            list(SeqIO.parse("GFA/corrupt_tag_name.gfa", "gfa1"))

    def test_corrupt_tag_type(self):
        """Check a GFA file with an incorrect tag type."""
        with self.assertWarns(BiopythonWarning):
            list(SeqIO.parse("GFA/corrupt_tag_type.gfa", "gfa1"))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
