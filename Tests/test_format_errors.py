# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Pin the exceptions raised for unknown, unsupported or malformed format names.

Bio.SeqIO, Bio.AlignIO, Bio.Align and Bio.Phylo each check a format name in
their own way and word their own messages.  These tables record, through the
public API only, what each entry point raises today, so that a change to how
formats are looked up cannot alter an exception type or message without a test
noticing.

Some rows pin quirks rather than designs: format(record, "txt") raises
KeyError, Bio.AlignIO accepts "123" where Bio.SeqIO rejects it, and
Bio.Align.write names an unwritable format in the case it was given.  Changing
any of them should be a deliberate decision, made together with its row here.
"""

import re
import unittest
from io import StringIO

import support

from Bio import Align
from Bio import AlignIO
from Bio import Phylo
from Bio import SeqIO
from Bio import StreamModeError
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

FASTA = support.DATA / "Fasta" / "f002"
SWISS = support.DATA / "SwissProt" / "F2CXE6.txt"
EMBOSS = support.DATA / "Emboss" / "needle.txt"
CLUSTAL = support.DATA / "Clustalw" / "opuntia.aln"
NEWICK = support.DATA / "Nexus" / "int_node_labels.nwk"

NOT_A_STRING = "Need a string for the file format (lower case)"
NO_FORMAT = "Format required (lower case string)"
# Python 3.14 reworded this message ("cannot use 'tuple' as a dict key
# (unhashable type: 'list')"), so only the common part is pinned.
UNHASHABLE = re.compile(re.escape("unhashable type: 'list'"))


class FormatErrorTestCase(unittest.TestCase):
    """Shared helper for checking a table of (format, exception, message)."""

    def check_table(self, call, table):
        """Check call(fmt) raises exactly the exception and message in each row."""
        for fmt, exception, message in table:
            with self.subTest(fmt=fmt):
                with self.assertRaises(Exception) as cm:
                    call(fmt)
                self.assertIs(type(cm.exception), exception)
                if isinstance(message, re.Pattern):
                    self.assertEqual(len(cm.exception.args), 1)
                    self.assertRegex(cm.exception.args[0], message)
                else:
                    self.assertEqual(cm.exception.args, (message,))


class SeqIOFormatErrors(FormatErrorTestCase):
    """Bio.SeqIO's entry points."""

    def test_parse(self):
        self.check_table(
            lambda fmt: SeqIO.parse(FASTA, fmt),
            [
                (None, TypeError, NOT_A_STRING),
                (42, TypeError, NOT_A_STRING),
                ("", ValueError, NO_FORMAT),
                ("FASTA", ValueError, "Format string 'FASTA' should be lower case"),
                ("123", ValueError, "Format string '123' should be lower case"),
                ("nope", ValueError, "Unknown format 'nope'"),
                (["fasta"], TypeError, NOT_A_STRING),
            ],
        )

    def test_write(self):
        record = SeqRecord(Seq("ACGT"), id="x", description="")
        self.check_table(
            lambda fmt: SeqIO.write([record], StringIO(), fmt),
            [
                (None, TypeError, NOT_A_STRING),
                (42, TypeError, NOT_A_STRING),
                ("", ValueError, NO_FORMAT),
                ("FASTA", ValueError, "Format string 'FASTA' should be lower case"),
                ("123", ValueError, "Format string '123' should be lower case"),
                ("nope", ValueError, "Unknown format 'nope'"),
                (["fasta"], TypeError, NOT_A_STRING),
                (
                    "swiss",
                    ValueError,
                    "Reading format 'swiss' is supported, but not writing",
                ),
                (
                    "emboss",
                    ValueError,
                    "Reading format 'emboss' is supported, but not writing",
                ),
            ],
        )

    def test_index(self):
        self.check_table(
            lambda fmt: SeqIO.index(str(FASTA), fmt),
            [
                (None, TypeError, NOT_A_STRING),
                (42, TypeError, NOT_A_STRING),
                ("", ValueError, NO_FORMAT),
                ("FASTA", ValueError, "Format string 'FASTA' should be lower case"),
                ("123", ValueError, "Format string '123' should be lower case"),
                ("nope", ValueError, "Unsupported format 'nope'"),
                (["fasta"], TypeError, NOT_A_STRING),
            ],
        )

    def test_index_db(self):
        missing = "Filenames to index and format required to build ':memory:'"
        self.check_table(
            lambda fmt: SeqIO.index_db(":memory:", str(FASTA), fmt),
            [
                (None, ValueError, missing),
                (42, TypeError, NOT_A_STRING),
                ("", ValueError, missing),
                ("FASTA", ValueError, "Format string 'FASTA' should be lower case"),
                ("123", ValueError, "Format string '123' should be lower case"),
                ("nope", ValueError, "Unsupported format 'nope'"),
                (["fasta"], TypeError, NOT_A_STRING),
            ],
        )

    def test_convert_input_format(self):
        self.check_table(
            lambda fmt: SeqIO.convert(FASTA, fmt, StringIO(), "fasta"),
            [
                (None, TypeError, NOT_A_STRING),
                (42, TypeError, NOT_A_STRING),
                ("", ValueError, NO_FORMAT),
                ("FASTA", ValueError, "Format string 'FASTA' should be lower case"),
                ("123", ValueError, "Format string '123' should be lower case"),
                ("nope", ValueError, "Unknown format 'nope'"),
                (["fasta"], TypeError, UNHASHABLE),
            ],
        )

    def test_convert_output_format(self):
        self.check_table(
            lambda fmt: SeqIO.convert(FASTA, "fasta", StringIO(), fmt),
            [
                (None, TypeError, NOT_A_STRING),
                (42, TypeError, NOT_A_STRING),
                ("", ValueError, NO_FORMAT),
                ("FASTA", ValueError, "Format string 'FASTA' should be lower case"),
                ("123", ValueError, "Format string '123' should be lower case"),
                ("nope", ValueError, "Unknown format 'nope'"),
                (["fasta"], TypeError, UNHASHABLE),
            ],
        )

    def test_seqrecord_format(self):
        record = SeqRecord(Seq("ACGT"), id="x", description="")
        binary = "Binary format sff cannot be used with SeqRecord format method"
        table = [
            ("txt", KeyError, "txt"),
            ("FASTA", KeyError, "FASTA"),
            ("sff", ValueError, binary),
        ]
        with self.subTest(call="format(record, fmt)"):
            self.check_table(lambda fmt: format(record, fmt), table)
        with self.subTest(call="record.format(fmt)"):
            self.check_table(record.format, table)
        with self.subTest(call="f-string"):
            self.check_table(lambda fmt: f"{record:{fmt}}", table)


class AlignIOFormatErrors(FormatErrorTestCase):
    """Bio.AlignIO's entry points, which consult Bio.SeqIO's tables too."""

    def test_parse(self):
        self.check_table(
            lambda fmt: next(AlignIO.parse(FASTA, fmt)),
            [
                ("FASTA", ValueError, "Format string 'FASTA' should be lower case"),
                ("123", ValueError, "Unknown format '123'"),
                ("nope", ValueError, "Unknown format 'nope'"),
            ],
        )

    def test_parse_through_either_table(self):
        """A SeqIO-only and an AlignIO-only format both parse."""
        lengths = [len(alignment) for alignment in AlignIO.parse(SWISS, "swiss")]
        self.assertEqual(lengths, [1])
        lengths = [len(alignment) for alignment in AlignIO.parse(EMBOSS, "emboss")]
        self.assertEqual(lengths, [2, 2, 2, 2, 2])

    def test_write(self):
        alignment = AlignIO.read(CLUSTAL, "clustal")
        self.check_table(
            lambda fmt: AlignIO.write([alignment], StringIO(), fmt),
            [
                ("FASTA", ValueError, "Format string 'FASTA' should be lower case"),
                ("123", ValueError, "Unknown format '123'"),
                ("nope", ValueError, "Unknown format 'nope'"),
                (
                    "swiss",
                    ValueError,
                    "Reading format 'swiss' is supported, but not writing",
                ),
                (
                    "emboss",
                    ValueError,
                    "Reading format 'emboss' is supported, but not writing",
                ),
            ],
        )


class AlignFormatErrors(FormatErrorTestCase):
    """Bio.Align's entry points, which ignore the case of a format name."""

    no_lower = "'NoneType' object has no attribute 'lower'"

    def test_parse(self):
        self.check_table(
            lambda fmt: Align.parse(CLUSTAL, fmt),
            [
                ("nope", ValueError, "Unknown file format nope"),
                ("NOPE", ValueError, "Unknown file format nope"),
                ("", ValueError, "Unknown file format "),
                (None, AttributeError, self.no_lower),
                (
                    "phylip-sequential",
                    ValueError,
                    "Unknown file format phylip-sequential",
                ),
            ],
        )

    def test_read(self):
        self.check_table(
            lambda fmt: Align.read(CLUSTAL, fmt),
            [
                ("nope", ValueError, "Unknown file format nope"),
                (None, AttributeError, self.no_lower),
            ],
        )

    def test_write(self):
        alignment = Align.read(CLUSTAL, "clustal")
        self.check_table(
            lambda fmt: Align.write(alignment, StringIO(), fmt),
            [
                ("nope", ValueError, "Unknown file format nope"),
                ("NOPE", ValueError, "Unknown file format nope"),
                ("", ValueError, "Unknown file format "),
                (None, AttributeError, self.no_lower),
                (
                    "phylip-sequential",
                    ValueError,
                    "Unknown file format phylip-sequential",
                ),
                (
                    "hhr",
                    ValueError,
                    "File writing has not yet been implemented for the hhr format",
                ),
                (
                    "HHR",
                    ValueError,
                    "File writing has not yet been implemented for the HHR format",
                ),
                ("bigbed", StreamModeError, "File must be opened in binary mode."),
            ],
        )

    def test_alignment_format(self):
        alignment = Align.read(CLUSTAL, "clustal")
        no_writer = (
            "Formatting alignments has not yet been implemented for the hhr format"
        )
        table = [
            ("nope", ValueError, "Unknown file format nope"),
            ("NOPE", ValueError, "Unknown file format nope"),
            ("phylip-sequential", ValueError, "Unknown file format phylip-sequential"),
            ("hhr", ValueError, no_writer),
            ("bigbed", ValueError, "bigbed is a binary file format"),
        ]
        with self.subTest(call="alignment.format(fmt)"):
            self.check_table(
                alignment.format, [*table, (None, AttributeError, self.no_lower)]
            )
        with self.subTest(call="format(alignment, fmt)"):
            self.check_table(lambda fmt: format(alignment, fmt), table)


class PhyloFormatErrors(FormatErrorTestCase):
    """Bio.Phylo's entry points, which match a format name exactly."""

    table = [
        ("nope", KeyError, "nope"),
        ("NEWICK", KeyError, "NEWICK"),
    ]

    def test_parse(self):
        self.check_table(lambda fmt: next(Phylo.parse(NEWICK, fmt)), self.table)

    def test_read(self):
        self.check_table(lambda fmt: Phylo.read(NEWICK, fmt), self.table)

    def test_write(self):
        tree = Phylo.read(NEWICK, "newick")
        self.check_table(lambda fmt: Phylo.write(tree, StringIO(), fmt), self.table)

    def test_convert_input_format(self):
        self.check_table(
            lambda fmt: Phylo.convert(NEWICK, fmt, StringIO(), "newick"), self.table
        )

    def test_convert_output_format(self):
        self.check_table(
            lambda fmt: Phylo.convert(NEWICK, "newick", StringIO(), fmt), self.table
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
