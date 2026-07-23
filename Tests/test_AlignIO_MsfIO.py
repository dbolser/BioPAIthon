# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for the Bio.AlignIO.MsfIO parser (GCG MSF format).

The real world MSF files under ``Tests/msf/`` exercise the happy path.  These
tests concentrate on malformed input, because that is what a parser is judged
on in practice.  Every expectation below is derived from the format as
described in ``Bio/AlignIO/MsfIO.py`` and in the GCG documentation quoted
there, not from observing what the current implementation happens to print.
"""

import unittest
import warnings
from io import StringIO

from Bio import AlignIO
from Bio import BiopythonParserWarning

# A ten column protein alignment, spelled out in the pieces the builder below
# assembles.  Ten columns fits inside a single fifty column block.
DEFAULT_HEADER = "!!AA_MULTIPLE_ALIGNMENT 1.0"
DEFAULT_MSF_LINE = "   demo.msf  MSF: 10  Type: P  January 6, 2000 15:41  Check: 0 .."
DEFAULT_NAMES = (
    " Name: Alpha  Len: 10  Check: 1755  Weight: 1.00",
    " Name: Beta  Len: 10  Check: 1035  Weight: 1.00",
)
DEFAULT_BLOCKS = (("Alpha  MKVLAAGIVL", "Beta   MKVLAAG.VL"),)

# The two ungapped sequences the default document above encodes.  The parser
# is documented to turn both "." and "~" gap characters into "-".
ALPHA_SEQ = "MKVLAAGIVL"
BETA_SEQ = "MKVLAAG-VL"


def build_msf(
    header=DEFAULT_HEADER,
    msf_line=DEFAULT_MSF_LINE,
    names=DEFAULT_NAMES,
    blocks=DEFAULT_BLOCKS,
    trailer=(),
):
    """Assemble a GCG MSF document from its parts and return it as a string.

    A blank line is emitted after the MSF/Type/Check line, after the closing
    ``//`` of the name block, and after every sequence block, which is the
    layout the GCG tools produce.
    """
    lines = [header, "", msf_line, ""]
    lines.extend(names)
    lines.extend(["", "//", ""])
    for block in blocks:
        lines.extend(block)
        lines.append("")
    lines.extend(trailer)
    return "\n".join(lines) + "\n"


def parse_msf(data):
    """Return the list of alignments parsed from an MSF document string."""
    return list(AlignIO.parse(StringIO(data), "msf"))


class MsfValidInput(unittest.TestCase):
    """Parsing of well formed, if unusual, MSF documents."""

    def test_minimal_alignment(self):
        """A two sequence, one block alignment round-trips its ids and gaps."""
        alignments = parse_msf(build_msf())
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(alignment.get_alignment_length(), 10)
        self.assertEqual([r.id for r in alignment], ["Alpha", "Beta"])
        self.assertEqual(str(alignment[0].seq), ALPHA_SEQ)
        self.assertEqual(str(alignment[1].seq), BETA_SEQ)
        self.assertEqual([r.annotations["weight"] for r in alignment], [1.00, 1.00])

    def test_tcoffee_oo_suffix_on_names(self):
        """T-COFFEE writes 'oo' after each name; it is not part of the id."""
        names = (
            " Name: Alpha oo  Len: 10  Check: 1755  Weight: 1.00",
            " Name: Beta oo  Len: 10  Check: 1035  Weight: 1.00",
        )
        alignments = parse_msf(build_msf(names=names))
        self.assertEqual(len(alignments), 1)
        self.assertEqual([r.id for r in alignments[0]], ["Alpha", "Beta"])

    def test_coordinate_header_lines(self):
        """Optional coordinate lines before each block are accepted and skipped."""
        msf_line = "   demo.msf  MSF: 60  Type: P  Check: 0 .."
        names = (
            " Name: Alpha  Len: 60  Check: 1755  Weight: 1.00",
            " Name: Beta  Len: 60  Check: 1035  Weight: 1.00",
        )
        blocks = (
            (
                "             1                                              50",
                "Alpha  " + "MKVLAAGIVL " * 4 + "MKVLAAGIVL",
                "Beta   " + "MKVLAAGIVL " * 4 + "MKVLAAGIVL",
            ),
            (
                "            51        60",
                "Alpha  ACGTACGTAC",
                "Beta   ACGTACGTAC",
            ),
        )
        alignments = parse_msf(build_msf(msf_line=msf_line, names=names, blocks=blocks))
        self.assertEqual(len(alignments), 1)
        self.assertEqual(alignments[0].get_alignment_length(), 60)
        self.assertEqual(str(alignments[0][0].seq), ALPHA_SEQ * 5 + "ACGTACGTAC")

    def test_tcoffee_double_blank_line_between_blocks(self):
        """T-COFFEE separates blocks with two blank lines rather than one."""
        msf_line = "   demo.msf  MSF: 60  Type: P  Check: 0 .."
        names = (
            " Name: Alpha  Len: 60  Check: 1755  Weight: 1.00",
            " Name: Beta  Len: 60  Check: 1035  Weight: 1.00",
        )
        blocks = (
            (
                "Alpha  " + "MKVLAAGIVL " * 4 + "MKVLAAGIVL",
                "Beta   " + "MKVLAAGIVL " * 4 + "MKVLAAGIVL",
                "",  # the builder adds a second blank line after this one
            ),
            (
                "Alpha  ACGTACGTAC",
                "Beta   ACGTACGTAC",
            ),
        )
        alignments = parse_msf(build_msf(msf_line=msf_line, names=names, blocks=blocks))
        self.assertEqual(len(alignments), 1)
        self.assertEqual(alignments[0].get_alignment_length(), 60)

    def test_short_sequence_omitted_from_final_block(self):
        """A sequence shorter than the alignment may simply stop, and is gap padded."""
        msf_line = "   demo.msf  MSF: 60  Type: P  Check: 0 .."
        names = (
            " Name: Alpha  Len: 60  Check: 1755  Weight: 1.00",
            " Name: Beta  Len: 10  Check: 1035  Weight: 1.00",
        )
        blocks = (
            (
                "Alpha  " + "MKVLAAGIVL " * 4 + "MKVLAAGIVL",
                "Beta   MKVLAAGIVL",
            ),
            (
                "Alpha  ACGTACGTAC",
                "",  # nothing at all for Beta in this block
            ),
        )
        data = build_msf(msf_line=msf_line, names=names, blocks=blocks)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            alignments = parse_msf(data)
        self.assertEqual(len(alignments), 1)
        self.assertEqual(alignments[0].get_alignment_length(), 60)
        self.assertEqual(str(alignments[0][1].seq), "MKVLAAGIVL" + "-" * 50)
        self.assertEqual(
            [w.category for w in caught if w.category is BiopythonParserWarning],
            [BiopythonParserWarning],
        )

    def test_two_concatenated_alignments(self):
        """Two MSF documents in one file give two alignments."""
        # The extra trailing blank lines between the two documents must be
        # skipped over rather than reported as unexpected content.
        alignments = parse_msf(build_msf(trailer=("", "")) + build_msf())
        self.assertEqual(len(alignments), 2)
        for alignment in alignments:
            self.assertEqual([r.id for r in alignment], ["Alpha", "Beta"])
            self.assertEqual(str(alignment[0].seq), ALPHA_SEQ)


class MsfBadHeader(unittest.TestCase):
    """Rejection of malformed MSF headers."""

    def test_unknown_first_line(self):
        """Only the three whitelisted header keywords may start the file."""
        data = build_msf(header="!!XX_MULTIPLE_ALIGNMENT 1.0")
        with self.assertRaises(ValueError) as cm:
            parse_msf(data)
        self.assertIn("not a known GCG MSF header", str(cm.exception))

    def test_missing_msf_line(self):
        """A file with no MSF/Type/Check line at all is rejected."""
        data = "!!AA_MULTIPLE_ALIGNMENT 1.0\n\nPileUp of: @nothing\n\n"
        with self.assertRaises(ValueError) as cm:
            parse_msf(data)
        self.assertIn("Reached end of file", str(cm.exception))

    def test_malformed_msf_line(self):
        """The MSF line must read 'MSF: <int> Type: <letter> ... Check: <int> ..'."""
        for msf_line in (
            # "Type:" keyword missing
            "   demo.msf  MSF: 10  Kind: P  Check: 0 ..",
            # neither "Check:" nor "CompCheck:"
            "   demo.msf  MSF: 10  Type: P  Checksum: 0 ..",
            # missing the trailing double dot
            "   demo.msf  MSF: 10  Type: P  Check: 0",
        ):
            with self.subTest(msf_line=msf_line):
                with self.assertRaises(ValueError) as cm:
                    parse_msf(build_msf(msf_line=msf_line))
                self.assertIn("GCG MSF header line should be", str(cm.exception))

    def test_compcheck_accepted(self):
        """EMBOSS writes 'CompCheck:' instead of 'Check:'; that is allowed."""
        msf_line = "   demo.msf  MSF: 10  Type: P  01/08/19 CompCheck: 8543 .."
        alignments = parse_msf(build_msf(msf_line=msf_line))
        self.assertEqual(len(alignments), 1)

    def test_bad_column_count(self):
        """The MSF value is the column count and must be a non-negative integer."""
        for msf_line in (
            "   demo.msf  MSF: ten  Type: P  Check: 0 ..",
            "   demo.msf  MSF: -5  Type: P  Check: 0 ..",
        ):
            with self.subTest(msf_line=msf_line):
                with self.assertRaises(ValueError) as cm:
                    parse_msf(build_msf(msf_line=msf_line))
                self.assertIn("for column count", str(cm.exception))

    def test_column_count_disagrees_with_name_lines(self):
        """The MSF column count must equal the longest Len: in the name block."""
        names = (
            " Name: Alpha  Len: 9  Check: 1755  Weight: 1.00",
            " Name: Beta  Len: 9  Check: 1035  Weight: 1.00",
        )
        with self.assertRaises(ValueError) as cm:
            parse_msf(build_msf(names=names))
        self.assertIn(
            "header said alignment length 10, but 2 of 2 sequences said Len: 9",
            str(cm.exception),
        )

    def test_bad_sequence_type(self):
        """Type must be P (protein) or N (nucleotide)."""
        msf_line = "   demo.msf  MSF: 10  Type: X  Check: 0 .."
        with self.assertRaises(ValueError) as cm:
            parse_msf(build_msf(msf_line=msf_line))
        self.assertIn("Type: X", str(cm.exception))


class MsfBadNameBlock(unittest.TestCase):
    """Rejection of malformed Name: lines and name blocks."""

    def test_duplicate_id(self):
        """Two sequences may not share an id."""
        names = (
            " Name: Alpha  Len: 10  Check: 1755  Weight: 1.00",
            " Name: Alpha  Len: 10  Check: 1035  Weight: 1.00",
        )
        with self.assertRaises(ValueError) as cm:
            parse_msf(build_msf(names=names))
        self.assertIn("Duplicated ID", str(cm.exception))

    def test_space_in_id(self):
        """Ids containing a space cannot be told apart from the sequence data."""
        names = (
            " Name: Alpha extra  Len: 10  Check: 1755  Weight: 1.00",
            " Name: Beta  Len: 10  Check: 1035  Weight: 1.00",
        )
        with self.assertRaises(NotImplementedError) as cm:
            parse_msf(build_msf(names=names))
        self.assertIn("Space in ID", str(cm.exception))

    def test_name_line_missing_fields(self):
        """A Name: line without Len:, Check: and Weight: is malformed."""
        names = (
            " Name: Alpha  Len: 10  Weight: 1.00",
            " Name: Beta  Len: 10  Check: 1035  Weight: 1.00",
        )
        with self.assertRaises(ValueError) as cm:
            parse_msf(build_msf(names=names))
        self.assertIn("Malformed GCG MSF name line", str(cm.exception))

    def test_truncated_before_end_of_name_block(self):
        """The name block must be closed by a // line before end of file."""
        data = (
            "!!AA_MULTIPLE_ALIGNMENT 1.0\n"
            "\n"
            "   demo.msf  MSF: 10  Type: P  Check: 0 ..\n"
            "\n"
            " Name: Alpha  Len: 10  Check: 1755  Weight: 1.00\n"
            " Name: Beta  Len: 10  Check: 1035  Weight: 1.00\n"
        )
        with self.assertRaises(ValueError) as cm:
            parse_msf(data)
        self.assertIn("End of file while looking for end of header", str(cm.exception))

    def test_truncated_after_slashes(self):
        """A file ending straight after // has no sequences to read."""
        data = (
            "!!AA_MULTIPLE_ALIGNMENT 1.0\n"
            "\n"
            "   demo.msf  MSF: 10  Type: P  Check: 0 ..\n"
            "\n"
            " Name: Alpha  Len: 10  Check: 1755  Weight: 1.00\n"
            " Name: Beta  Len: 10  Check: 1035  Weight: 1.00\n"
            "\n"
            "//\n"
        )
        with self.assertRaises(ValueError) as cm:
            parse_msf(data)
        self.assertIn("End of file after // line", str(cm.exception))

    def test_no_blank_line_after_slashes(self):
        """The // line must be followed by a blank line, not by sequence data."""
        data = (
            "!!AA_MULTIPLE_ALIGNMENT 1.0\n"
            "\n"
            "   demo.msf  MSF: 10  Type: P  Check: 0 ..\n"
            "\n"
            " Name: Alpha  Len: 10  Check: 1755  Weight: 1.00\n"
            " Name: Beta  Len: 10  Check: 1035  Weight: 1.00\n"
            "\n"
            "//\n"
            "Alpha  MKVLAAGIVL\n"
            "Beta   MKVLAAGIVL\n"
        )
        with self.assertRaises(ValueError) as cm:
            parse_msf(data)
        self.assertIn("expected blank line before sequences", str(cm.exception))


class MsfBadSequenceBlocks(unittest.TestCase):
    """Rejection of malformed sequence blocks."""

    MSF_LINE_60 = "   demo.msf  MSF: 60  Type: P  Check: 0 .."
    NAMES_60 = (
        " Name: Alpha  Len: 60  Check: 1755  Weight: 1.00",
        " Name: Beta  Len: 60  Check: 1035  Weight: 1.00",
    )
    FULL_BLOCK = (
        "Alpha  " + "MKVLAAGIVL " * 4 + "MKVLAAGIVL",
        "Beta   " + "MKVLAAGIVL " * 4 + "MKVLAAGIVL",
    )

    def build60(self, blocks):
        """Return a sixty column document made of the given sequence blocks."""
        return build_msf(msf_line=self.MSF_LINE_60, names=self.NAMES_60, blocks=blocks)

    def test_truncated_before_last_block(self):
        """Running out of file part way through the blocks is an error."""
        with self.assertRaises(ValueError) as cm:
            parse_msf(self.build60((self.FULL_BLOCK,)))
        self.assertIn("End of file where expecting sequence data", str(cm.exception))

    def test_coordinate_line_wrong_start(self):
        """A coordinate line must start at the next uncompleted column."""
        blocks = (("             5        50",) + self.FULL_BLOCK,)
        with self.assertRaises(ValueError) as cm:
            parse_msf(self.build60(blocks))
        self.assertIn("coordinate line starting 1", str(cm.exception))

    def test_coordinate_line_non_numeric_start(self):
        """A first token that is neither the first id nor a number is an error."""
        blocks = (("             one        50",) + self.FULL_BLOCK,)
        with self.assertRaises(ValueError) as cm:
            parse_msf(self.build60(blocks))
        self.assertIn("coordinate line starting 1", str(cm.exception))

    def test_coordinate_line_wrong_end(self):
        """The second coordinate must be the last column of this block."""
        blocks = (("             1        99",) + self.FULL_BLOCK,)
        with self.assertRaises(ValueError) as cm:
            parse_msf(self.build60(blocks))
        self.assertIn("coordinate line 1 to 50", str(cm.exception))

    def test_coordinate_line_non_numeric_end(self):
        """A non-numeric end coordinate is an error."""
        blocks = (("             1        fifty",) + self.FULL_BLOCK,)
        with self.assertRaises(ValueError) as cm:
            parse_msf(self.build60(blocks))
        self.assertIn("coordinate line 1 to 50", str(cm.exception))

    def test_coordinate_line_too_many_fields(self):
        """A coordinate line carries at most a start and an end coordinate."""
        blocks = (("             1     25     50",) + self.FULL_BLOCK,)
        with self.assertRaises(ValueError) as cm:
            parse_msf(self.build60(blocks))
        self.assertIn("coordinate line 1 to 50", str(cm.exception))

    def test_coordinate_line_too_many_fields_final_block(self):
        """A three field coordinate line is wrong even if the start matches.

        On the final block of a 51 column alignment the start coordinate, 51,
        is also the end coordinate, so only the count of fields on the line
        can tell the parser that the line is malformed.
        """
        msf_line = "   demo.msf  MSF: 51  Type: P  Check: 0 .."
        names = (
            " Name: Alpha  Len: 51  Check: 1755  Weight: 1.00",
            " Name: Beta  Len: 51  Check: 1035  Weight: 1.00",
        )
        blocks = (
            self.FULL_BLOCK,
            ("            51    51    51", "Alpha  M", "Beta   M"),
        )
        with self.assertRaises(ValueError) as cm:
            parse_msf(build_msf(msf_line=msf_line, names=names, blocks=blocks))
        self.assertIn("coordinate line 51 to 51", str(cm.exception))

    def test_missing_sequence_line(self):
        """A blank line in place of a full length sequence is an error."""
        blocks = ((self.FULL_BLOCK[0], ""),)
        with self.assertRaises(ValueError) as cm:
            parse_msf(self.build60(blocks))
        self.assertIn("Expected sequence for Beta", str(cm.exception))

    def test_wrong_name_in_block(self):
        """Sequence lines must appear in the order given by the name block."""
        blocks = ((self.FULL_BLOCK[0], "Gamma  " + "MKVLAAGIVL " * 4 + "MKVLAAGIVL"),)
        with self.assertRaises(ValueError) as cm:
            parse_msf(self.build60(blocks))
        self.assertIn("Expected sequence for 'Beta'", str(cm.exception))

    def test_missing_blank_line_between_blocks(self):
        """Sequence blocks must be separated by a blank line."""
        data = self.build60(
            (self.FULL_BLOCK + ("Alpha  ACGTACGTAC", "Beta   ACGTACGTAC"),)
        )
        with self.assertRaises(ValueError) as cm:
            parse_msf(data)
        self.assertIn("Expected blank line", str(cm.exception))

    def test_junk_after_alignment(self):
        """Anything other than blank lines or a new header after an alignment."""
        with self.assertRaises(ValueError) as cm:
            parse_msf(build_msf(trailer=("not an MSF header",)))
        self.assertIn("Unexpected line after GCG MSF alignment", str(cm.exception))

    def test_sequences_longer_than_declared(self):
        """Sequence data exceeding the declared column count is an error."""
        # Both the MSF line and the Name: lines say ten columns, but each
        # sequence line carries twenty residues.
        blocks = (("Alpha  MKVLAAGIVL MKVLAAGIVL", "Beta   MKVLAAGIVL MKVLAAGIVL"),)
        with self.assertRaises(ValueError) as cm:
            parse_msf(build_msf(blocks=blocks))
        self.assertIn(
            "headers said alignment length 10, but have 20", str(cm.exception)
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
