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
