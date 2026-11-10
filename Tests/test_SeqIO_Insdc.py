# Copyright 2013 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for SeqIO Insdc module."""

import unittest
import warnings
from io import StringIO

import support
from seq_tests_common import SeqRecordTestBaseClass
from test_SeqIO import SeqIOConverterTestBaseClass

from Bio import BiopythonParserWarning
from Bio import BiopythonWarning
from Bio import GenBank
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import SimpleLocation
from Bio.SeqRecord import SeqRecord


class TestEmbl(unittest.TestCase):
    def test_annotation1(self):
        """Check parsing of annotation from EMBL files (1)."""
        record = SeqIO.read(support.DATA / "EMBL" / "TRBG361.embl", "embl")
        self.assertEqual(len(record), 1859)
        # Single keyword:
        self.assertEqual(record.annotations["keywords"], ["beta-glucosidase"])
        self.assertEqual(record.annotations["topology"], "linear")

    def test_annotation2(self):
        """Check parsing of annotation from EMBL files (2)."""
        record = SeqIO.read(support.DATA / "EMBL" / "DD231055_edited.embl", "embl")
        self.assertEqual(len(record), 315)
        # Multiple keywords:
        self.assertEqual(
            record.annotations["keywords"],
            [
                "JP 2005522996-A/12",
                "test-data",
                "lot and lots of keywords for this example",
                "multi-line keywords",
            ],
        )
        self.assertEqual(record.annotations["topology"], "linear")

    def test_annotation3(self):
        """Check parsing of annotation from EMBL files (3)."""
        record = SeqIO.read(support.DATA / "EMBL" / "AE017046.embl", "embl")
        self.assertEqual(len(record), 9609)
        # TODO: Should this be an empty list, or simply absent?
        self.assertEqual(record.annotations["keywords"], [""])
        self.assertEqual(record.annotations["topology"], "circular")

    def test_annotation4(self):
        """Check parsing of annotation from EMBL files (4)."""
        with self.assertWarns(BiopythonParserWarning):
            record = SeqIO.read(support.DATA / "EMBL" / "location_wrap.embl", "embl")
        self.assertEqual(len(record), 120)
        self.assertNotIn("keywords", record.annotations)
        # The ID line has the topology as unspecified:
        self.assertNotIn("topology", record.annotations)

    def test_writing_empty_qualifiers(self):
        f = SeqFeature(
            SimpleLocation(5, 20, strand=+1),
            type="region",
            qualifiers={"empty": None, "zero": 0, "one": 1, "text": "blah"},
        )
        record = SeqRecord(Seq("A" * 100), "dummy", features=[f])
        record.annotations["molecule_type"] = "DNA"
        gbk = record.format("gb")
        self.assertIn(" /empty\n", gbk)
        self.assertIn(" /zero=0\n", gbk)
        self.assertIn(" /one=1\n", gbk)
        self.assertIn(' /text="blah"\n', gbk)

    def test_warn_on_writing_nonstandard_feature_key(self):
        f = SeqFeature(
            SimpleLocation(5, 20, strand=+1),
            type="a" * 16,
            qualifiers={"empty": None, "zero": 0, "one": 1, "text": "blah"},
        )
        record = SeqRecord(Seq("A" * 100), "dummy", features=[f])
        record.annotations["molecule_type"] = "DNA"
        with self.assertWarns(BiopythonWarning):
            record.format("gb")

    def test_warn_on_writing_nonstandard_qualifier_key(self):
        f = SeqFeature(
            SimpleLocation(5, 20, strand=+1),
            type="region",
            qualifiers={"a" * 21: "test"},
        )
        record = SeqRecord(Seq("A" * 100), "dummy", features=[f])
        record.annotations["molecule_type"] = "DNA"
        with self.assertWarns(BiopythonWarning):
            record.format("gb")


class TestValuelessQualifiers(unittest.TestCase):
    """Qualifiers that take no value, like /pseudo, are written bare.

    The parsers store a bare /pseudo as [""], which the writers used to
    write as /pseudo="", invalid INSDC (upstream issues #902 and #4205).
    """

    indent = {"gb": " " * 21, "embl": "FT" + " " * 19, "imgt": "FT" + " " * 23}

    def check_round_trip(self, record, key, count):
        """Write record in each format; check key is bare and survives."""
        for fmt, indent in self.indent.items():
            with self.subTest(fmt=fmt):
                text = record.format(fmt)
                lines = text.splitlines()
                self.assertEqual(lines.count(f"{indent}/{key}"), count)
                self.assertNotIn(f"/{key}=", text)
                new = SeqIO.read(StringIO(text), fmt)
                self.assertEqual(
                    [f.qualifiers for f in new.features],
                    [f.qualifiers for f in record.features],
                )

    def test_pseudo(self):
        record = SeqIO.read(support.DATA / "GenBank" / "blank_seq.gb", "gb")
        self.assertEqual(record.features[3].qualifiers["pseudo"], [""])
        self.check_round_trip(record, "pseudo", 1)

    def test_trans_splicing(self):
        record = SeqIO.read(support.DATA / "GenBank" / "NC_000932.gb", "gb")
        self.check_round_trip(record, "trans_splicing", 4)

    def test_empty_string(self):
        f = SeqFeature(
            SimpleLocation(5, 20, strand=+1),
            type="CDS",
            qualifiers={"pseudo": "", "germline": [""]},
        )
        record = SeqRecord(Seq("A" * 100), "dummy", features=[f])
        record.annotations["molecule_type"] = "DNA"
        for fmt, indent in self.indent.items():
            with self.subTest(fmt=fmt):
                lines = record.format(fmt).splitlines()
                self.assertIn(f"{indent}/pseudo", lines)
                self.assertIn(f"{indent}/germline", lines)

    def test_empty_note_unchanged(self):
        """An empty value for a qualifier that takes one is still quoted."""
        f = SeqFeature(
            SimpleLocation(5, 20, strand=+1),
            type="misc_feature",
            qualifiers={"note": [""]},
        )
        record = SeqRecord(Seq("A" * 100), "dummy", features=[f])
        record.annotations["molecule_type"] = "DNA"
        for fmt, indent in self.indent.items():
            with self.subTest(fmt=fmt):
                lines = record.format(fmt).splitlines()
                self.assertIn(f'{indent}/note=""', lines)
                self.assertNotIn(f"{indent}/note", lines)


class TestNoSpaceQualifiers(unittest.TestCase):
    """Wrapped /anticodon and /transl_except values are joined with no space.

    The NCBI breaks these unquoted values after a comma and said that a
    parser must append the next line with no space (upstream issue #2112).
    The parsers used to join the lines with a space, and the writers broke
    a long value before its closing bracket, which then read back with a
    space too.
    """

    anticodon = "(pos:complement(6180084..6180086),aa:Ile,seq:gat)"
    transl_except = "(pos:complement(1123552..1123554),aa:Sec)"

    # The /anticodon lines are copied from NC_018080. The /transl_except wrap
    # is made up: the NCBI said a break anywhere is legal, but the one
    # /transl_except in each of NC_018080 and NC_008463 fits on one line.
    genbank = """\
LOCUS       TEST                     100 bp    DNA              UNK 01-JAN-2000
DEFINITION  Wrapped qualifiers.
ACCESSION   TEST
VERSION     TEST.1
FEATURES             Location/Qualifiers
     source          1..100
     tRNA            complement(10..80)
                     /anticodon=(pos:complement(6180084..6180086),aa:Ile,
                     seq:gat)
                     /note="a wrapped
                     note"
     CDS             complement(1..99)
                     /transl_except=(pos:complement(1123552..1123554),
                     aa:Sec)
                     /protein_id="XP_000001.1"
                     /translation="MU"
ORIGIN
        1 aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa
       61 aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa
//
"""

    embl = """\
ID   TEST; SV 1; linear; DNA; STD; UNC; 100 BP.
XX
AC   TEST;
XX
DE   Wrapped qualifiers
XX
FH   Key             Location/Qualifiers
FH
FT   source          1..100
FT   tRNA            complement(10..80)
FT                   /anticodon=(pos:complement(6180084..6180086),aa:Ile,
FT                   seq:gat)
FT                   /note="a wrapped
FT                   note"
FT   CDS             complement(1..99)
FT                   /transl_except=(pos:complement(1123552..1123554),
FT                   aa:Sec)
FT                   /protein_id="XP_000001.1"
FT                   /translation="MU"
XX
SQ   Sequence 100 BP; 100 A; 0 C; 0 G; 0 T; 0 other;
     aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa        60
     aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa                             100
//
"""

    def test_read(self):
        for fmt in ("genbank", "embl"):
            with self.subTest(fmt=fmt):
                record = SeqIO.read(StringIO(getattr(self, fmt)), fmt)
                trna, cds = record.features[1:]
                self.assertEqual(trna.qualifiers["anticodon"], [self.anticodon])
                self.assertEqual(cds.qualifiers["transl_except"], [self.transl_except])
                # Other wrapped values are still joined with a space:
                self.assertEqual(trna.qualifiers["note"], ["a wrapped note"])

    def test_read_cds(self):
        for fmt in ("genbank", "embl"):
            with self.subTest(fmt=fmt):
                text = getattr(self, fmt)
                record = SeqIO.read(StringIO(text), f"{fmt}-cds")
                self.assertEqual(
                    record.annotations["transl_except"], self.transl_except
                )

    def test_read_quoted(self):
        """A quoted wrapped value, as BioPerl once wrote, is joined the same way."""
        text = self.genbank.replace("/anticodon=(", '/anticodon="(')
        text = text.replace("seq:gat)", 'seq:gat)"')
        record = SeqIO.read(StringIO(text), "genbank")
        self.assertEqual(record.features[1].qualifiers["anticodon"], [self.anticodon])

    def test_read_old_writer_output(self):
        """A value the old writers broke inside a number reads back whole."""
        indent = " " * 21
        # As Biopython 1.88 wrote it, broken at column 80:
        text = self.genbank.replace(
            f"{indent}/transl_except=(pos:complement(1123552..1123554),\n{indent}aa:Sec)\n",
            f"{indent}/transl_except=(pos:complement(join(1234567..1234568,123457\n"
            f"{indent}0)),aa:Sec)\n",
        )
        record = SeqIO.read(StringIO(text), "genbank")
        self.assertEqual(
            record.features[2].qualifiers["transl_except"],
            ["(pos:complement(join(1234567..1234568,1234570)),aa:Sec)"],
        )

    def test_write(self):
        """Long values break after a comma and read back unchanged."""
        anticodon = "(pos:complement(1234567..1234569),aa:Leu,seq:tag)"
        transl_except = "(pos:complement(join(1234567..1234568,1234570)),aa:Sec)"
        record = SeqRecord(Seq("A" * 100), "dummy")
        record.annotations["molecule_type"] = "DNA"
        record.features = [
            SeqFeature(
                SimpleLocation(9, 80, strand=-1),
                type="tRNA",
                qualifiers={"anticodon": [anticodon]},
            ),
            SeqFeature(
                SimpleLocation(0, 99, strand=-1),
                type="CDS",
                qualifiers={"transl_except": [transl_except]},
            ),
        ]
        indents = {"gb": " " * 21, "embl": "FT" + " " * 19, "imgt": "FT" + " " * 23}
        for fmt, indent in indents.items():
            with self.subTest(fmt=fmt):
                text = record.format(fmt)
                lines = text.splitlines()
                self.assertLessEqual(max(len(line) for line in lines), 80)
                self.assertIn(
                    f"{indent}/anticodon=(pos:complement(1234567..1234569),aa:Leu,",
                    lines,
                )
                self.assertIn(f"{indent}seq:tag)", lines)
                new = SeqIO.read(StringIO(text), fmt)
                self.assertEqual(new.features[0].qualifiers["anticodon"], [anticodon])
                self.assertEqual(
                    new.features[1].qualifiers["transl_except"], [transl_except]
                )

    def test_write_comma_at_width_limit(self):
        """A comma up to column 80 can end the line; one in column 81 cannot."""
        indent = " " * 21
        start = f"{indent}/transl_except=(pos:join("  # 46 columns
        for column in (79, 80, 81):
            n = column - 50  # puts the first comma in this column
            value = f"(pos:join({'1' * n}..2,3..4),aa:Sec)"
            with self.subTest(column=column):
                record = SeqRecord(Seq("A" * 100), "dummy")
                record.annotations["molecule_type"] = "DNA"
                record.features = [
                    SeqFeature(
                        SimpleLocation(0, 99),
                        type="CDS",
                        qualifiers={"transl_except": [value]},
                    )
                ]
                text = record.format("gb")
                lines = text.splitlines()
                if column <= 80:
                    expected = [f"{start}{'1' * n}..2,", f"{indent}3..4),aa:Sec)"]
                else:  # no comma fits, so break at the width limit
                    expected = [f"{start}{'1' * n}..2", f"{indent},3..4),aa:Sec)"]
                i = next(i for i, line in enumerate(lines) if line.startswith(start))
                self.assertEqual(lines[i : i + 2], expected)
                self.assertEqual(len(lines[i]), min(column, 80))
                new = SeqIO.read(StringIO(text), "gb")
                self.assertEqual(new.features[0].qualifiers["transl_except"], [value])

    def test_record_parser_and_writer(self):
        """Bio.GenBank's record parser and its writer agree."""
        record = GenBank.read(StringIO(self.genbank))
        trna, cds = record.features[1:]
        self.assertEqual(trna.qualifiers[0].value, self.anticodon)
        self.assertEqual(cds.qualifiers[0].value, self.transl_except)
        lines = str(record).splitlines()
        # Same layout as the NCBI's, with the comma ending the first line:
        self.assertIn(
            "                     /anticodon=(pos:complement(6180084..6180086),aa:Ile,",
            lines,
        )
        self.assertIn("                     seq:gat)", lines)
        self.assertLessEqual(max(len(line) for line in lines), 80)
        new = GenBank.read(StringIO(str(record)))
        self.assertEqual(new.features[1].qualifiers[0].value, self.anticodon)
        self.assertEqual(new.features[2].qualifiers[0].value, self.transl_except)


class TestInsertionBetweenBases(unittest.TestCase):
    """Ensembl's 7..6, an insertion between bases 6 and 7, is read as 6^7.

    Ensembl's GenBank and EMBL exports write the insertion this way with the
    sequence's real length. The parsers took it for a feature spanning the
    origin of a sequence that is not circular, so they warned and set the
    location to None, and writing the record then failed (upstream issue
    #4969). On a circular sequence 7..6 still spans the origin.
    """

    genbank = """\
LOCUS       TEST                      10 bp    DNA     {topology:8} UNK 01-JAN-1980
DEFINITION  Insertion.
ACCESSION   TEST
VERSION     TEST
FEATURES             Location/Qualifiers
     source          1..10
     variation       7..6
                     /replace="-/T"
ORIGIN
        1 aaagaaaggc
//
"""

    embl = """\
ID   TEST; SV 1; {topology}; DNA; STD; UNC; 10 BP.
XX
AC   TEST;
XX
DE   Insertion
XX
FH   Key             Location/Qualifiers
FH
FT   source          1..10
FT   variation       7..6
FT                   /replace="-/T"
XX
SQ   Sequence 10 BP; 6 A; 1 C; 3 G; 0 T; 0 other;
     aaagaaaggc                                                               10
//
"""

    def read(self, fmt, topology):
        text = getattr(self, fmt).format(topology=topology)
        return SeqIO.read(StringIO(text), fmt)

    def test_not_circular(self):
        """With no topology, or a linear one, the insertion is [6:6]."""
        for fmt, topology in [
            ("genbank", ""),
            ("genbank", "linear"),
            ("embl", "linear"),
        ]:
            with self.subTest(fmt=fmt, topology=topology):
                with warnings.catch_warnings():
                    warnings.simplefilter("error")
                    record = self.read(fmt, topology)
                location = record.features[1].location
                self.assertEqual(location, SimpleLocation(6, 6, strand=1))
                text = record.format(fmt)
                self.assertIn("variation       6^7\n", text)
                self.assertEqual(
                    SeqIO.read(StringIO(text), fmt).features, record.features
                )

    def test_circular(self):
        """On a circular sequence 7..6 still spans the origin."""
        for fmt in ["genbank", "embl"]:
            with self.subTest(fmt=fmt):
                with self.assertWarnsRegex(BiopythonParserWarning, "origin wrapping"):
                    record = self.read(fmt, "circular")
                self.assertEqual(
                    record.features[1].location,
                    SimpleLocation(6, 10, strand=1) + SimpleLocation(0, 6, strand=1),
                )
                text = record.format(fmt)
                self.assertIn("variation       join(7..10,1..6)\n", text)


class TestWriterErrors(unittest.TestCase):
    """Records the INSDC writers cannot represent raise ValueError."""

    def make_record(self, name="dummy", seq="ACGT"):
        record = SeqRecord(Seq(seq), id=name, name=name)
        record.annotations["molecule_type"] = "DNA"
        return record

    def check(self, record, msg, formats=("genbank", "embl")):
        for fmt in formats:
            with self.subTest(fmt=fmt):
                with self.assertRaises(ValueError) as cm:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", BiopythonWarning)
                        record.format(fmt)
                self.assertEqual(str(cm.exception), msg)

    def test_location_ref_db(self):
        record = self.make_record()
        location = SimpleLocation(0, 2, ref="X", ref_db="db")
        record.features.append(SeqFeature(location, type="misc_feature"))
        self.check(
            record,
            "Location db:X[0:2] has ref_db 'db', "
            "which INSDC location strings cannot represent",
        )

    def test_feature_without_type(self):
        record = self.make_record()
        record.features.append(SeqFeature(SimpleLocation(0, 2)))
        self.check(record, "Cannot write a feature with no type, at location [0:2]")

    def test_feature_without_location(self):
        # The parsers leave None, with a warning, for a location they cannot read
        record = self.make_record()
        record.features.append(SeqFeature(None, type="variation"))
        self.check(
            record,
            "Cannot write a 'variation' feature with no location "
            "(its .location is None)",
            formats=("genbank", "embl", "imgt"),
        )

    def test_annotation_list(self):
        record = self.make_record()
        record.annotations["organism"] = ["Homo sapiens", "Mus musculus"]
        self.check(
            record,
            "Expected a single value for annotation 'organism', "
            "not ['Homo sapiens', 'Mus musculus']",
        )

    def test_segment_list(self):
        record = self.make_record()
        record.annotations["segment"] = ["1 of 2", "2 of 2"]
        self.check(
            record,
            "Expected a single value for annotation 'segment', "
            "not ['1 of 2', '2 of 2']",
            formats=("genbank",),
        )

    def test_locus_with_whitespace(self):
        for name in ("dummy ", "dummy\n"):
            with self.subTest(name=name):
                record = self.make_record(name)
                with self.assertRaises(ValueError) as cm:
                    record.format("genbank")
                self.assertTrue(
                    str(cm.exception).startswith(
                        f"LOCUS line does not contain the locus {name!r} and "
                        "length 4 at the expected positions:\n"
                    ),
                    str(cm.exception),
                )

    def test_locus_and_length_too_long(self):
        # 16 character name and 12 digit length leave no space between them
        record = self.make_record("ABCDEFGHIJKLMNOP")
        record.seq = Seq(None, 10**11)
        self.check(
            record,
            "Locus name 'ABCDEFGHIJKLMNOP' and sequence length 100000000000 "
            "do not fit in the LOCUS line",
            formats=("genbank",),
        )


class TestEmblRewrite(SeqRecordTestBaseClass):
    def check_rewrite(self, filename):
        old = SeqIO.read(filename, "embl")

        # TODO - Check these properties:
        old.dbxrefs = []
        old.annotations["accessions"] = old.annotations["accessions"][:1]
        del old.annotations["references"]

        buffer = StringIO()
        self.assertEqual(1, SeqIO.write(old, buffer, "embl"))
        buffer.seek(0)
        new = SeqIO.read(buffer, "embl")

        self.compare_record(old, new)

    def test_annotation1(self):
        """Check writing-and-parsing EMBL file (1)."""
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            self.check_rewrite(support.DATA / "EMBL" / "TRBG361.embl")

    def test_annotation2(self):
        """Check writing-and-parsing EMBL file (2)."""
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            self.check_rewrite(support.DATA / "EMBL" / "DD231055_edited.embl")

    def test_annotation3(self):
        """Check writing-and-parsing EMBL file (3)."""
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            self.check_rewrite(support.DATA / "EMBL" / "AE017046.embl")


class ConvertTestsInsdc(SeqIOConverterTestBaseClass):
    def test_conversion(self):
        """Test format conversion by SeqIO.write/SeqIO.parse and SeqIO.convert."""
        tests = [
            (support.DATA / "EMBL" / "U87107.embl", "embl"),
            (support.DATA / "EMBL" / "TRBG361.embl", "embl"),
            (support.DATA / "GenBank" / "NC_005816.gb", "gb"),
            (support.DATA / "GenBank" / "cor6_6.gb", "genbank"),
        ]
        for filename, fmt in tests:
            for in_format, out_format in self.formats:
                if in_format != fmt:
                    continue
                self.check_conversion(filename, in_format, out_format)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
