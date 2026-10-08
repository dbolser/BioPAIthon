# Copyright 2015 by Kai Blin.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for EMBL module (using unittest framework)."""

import unittest
import warnings
from io import StringIO

import support

from Bio import BiopythonParserWarning
from Bio import SeqIO


class EMBLTests(unittest.TestCase):
    def test_embl_content_after_co(self):
        """Test a ValueError is thrown by content after a CO line."""

        def parse_content_after_co():
            rec = SeqIO.read(support.DATA / "EMBL" / "xx_after_co.embl", "embl")

        self.assertRaises(ValueError, parse_content_after_co)

        try:
            parse_content_after_co()
        except ValueError as e:
            self.assertEqual(str(e), "Unexpected content after SQ or CO line: 'XX'")
        else:
            self.assertTrue(
                False,
                "Error message without explanation raised by content after CO line",
            )

    def test_embl_0_line(self):
        """Test SQ line with 0 length sequence."""
        # Biopython 1.67 and older would parse this file with a warning:
        # 'Expected sequence length 1740, found 1744 (TIR43YW1_CE).' and
        # the coordinates 1740 added to the sequence as four extra letters.
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            rec = SeqIO.read(support.DATA / "EMBL" / "embl_with_0_line.embl", "embl")
            self.assertEqual(
                len(w),
                0,
                "Unexpected parser warnings: "
                + "\n".join(str(warn.message) for warn in w),
            )
            self.assertEqual(len(rec), 1740)

    def test_embl_no_coords(self):
        """Test sequence lines without coordinates."""
        # Biopython 1.68, 1.69 and 1.70 would ignore these lines
        # giving an unknown sequence!
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", BiopythonParserWarning)
            rec = SeqIO.read(support.DATA / "EMBL" / "101ma_no_coords.embl", "embl")
            self.assertTrue(w, "Expected parser warning")
            self.assertEqual(
                [str(_.message) for _ in w],
                ["EMBL sequence line missing coordinates"] * 3,
            )
            self.assertEqual(len(rec), 154)
            self.assertEqual(rec.seq[:10], "MVLSEGEWQL")
            self.assertEqual(rec.seq[-10:], "AKYKELGYQG")

    def test_embl_wrong_dr_line(self):
        """Test files with wrong DR lines."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", BiopythonParserWarning)
            record = SeqIO.read(support.DATA / "EMBL" / "RepBase23.02.embl", "embl")
            self.assertTrue(w, "Expected parser warning")
            self.assertEqual(
                [str(_.message) for _ in w], ["Malformed DR line in EMBL file."]
            )


class MalformedInputTests(unittest.TestCase):
    """Malformed EMBL and IMGT input raises ValueError.

    These were asserts, so under python -O the checks vanished; a truncated
    IMGT feature then made the parser loop forever.
    """

    def read_modified(self, filename, fmt, old, new):
        with open(support.DATA / "EMBL" / filename) as handle:
            data = handle.read()
        self.assertIn(old, data)
        return SeqIO.read(StringIO(data.replace(old, new, 1)), fmt)

    def test_length_unit(self):
        with self.assertRaisesRegex(ValueError, "Expected sequence length in BP"):
            self.read_modified("AAA03323.embl", "embl", "1545 BP.", "1545 XX.")

    def test_length_field(self):
        with self.assertRaisesRegex(ValueError, "Invalid sequence length string"):
            self.read_modified("AAA03323.embl", "embl", "1545 BP.", "1545.")

    def test_imgt_feature_line_not_ft(self):
        with self.assertRaisesRegex(ValueError, "Expected an FT line"):
            self.read_modified(
                "A04195.imgt", "imgt", "FT   INIT-CODON", "CC   INIT-CODON"
            )

    def test_imgt_blank_line_in_feature(self):
        with self.assertRaisesRegex(
            ValueError, "Expected an FT line in the 'L-REGION'"
        ):
            self.read_modified("A04195.imgt", "imgt", "/partial\n", "/partial\n\n")

    def test_imgt_truncated_in_feature(self):
        with open(support.DATA / "EMBL" / "A04195.imgt") as handle:
            data = handle.read()
        data = data[: data.index("/partial\n") + len("/partial\n")]
        with self.assertRaisesRegex(ValueError, "Premature end of file in the"):
            SeqIO.read(StringIO(data), "imgt")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
