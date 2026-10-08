# Copyright 2009-2013 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Additional unit tests for Bio.SeqIO.QualityIO (covering FASTQ and QUAL)."""

import array
import os
import pickle
import unittest
import warnings
from io import BytesIO
from io import StringIO

import support
from test_SeqIO import SeqIOConverterTestBaseClass
from test_SeqIO import SeqIOTestBaseClass

from Bio import BiopythonParserWarning
from Bio import BiopythonWarning
from Bio import SeqIO
from Bio.Data.IUPACData import ambiguous_dna_letters
from Bio.Data.IUPACData import ambiguous_rna_letters
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from Bio.Seq import UndefinedSequenceError
from Bio.SeqIO import QualityIO
from Bio.SeqRecord import SeqRecord


class QualityIOTestBaseClass(SeqIOTestBaseClass):
    def compare_record(self, old, new, fmt=None, msg=None):
        """Quality-aware SeqRecord comparison.

        This will check the mapping between Solexa and PHRED scores.
        It knows to ignore records with undefined sequences  for string
        matching (i.e. QUAL files) via the base class.
        """
        super().compare_record(old, new, msg=None)
        if fmt in ["fastq-solexa", "fastq-illumina"]:
            truncate = 62
        elif fmt in ["fastq", "fastq-sanger"]:
            truncate = 93
        else:
            assert fmt in ["fasta", "qual", "phd", "sff", "tab", None]
            truncate = None
        for keyword in ("phred_quality", "solexa_quality"):
            q_old = old.letter_annotations.get(keyword)
            q_new = new.letter_annotations.get(keyword)
            if q_old is None or q_new is None:
                continue
            if truncate is not None and q_old != q_new:
                q_old = [min(q, truncate) for q in q_old]
                q_new = [min(q, truncate) for q in q_new]
            err_msg = f"mismatch in {keyword}"
            if msg is not None:
                err_msg = f"{msg}: {err_msg}"
            self.assertEqual(q_old, q_new, msg=err_msg)

        q_old = old.letter_annotations.get("phred_quality")
        q_new = new.letter_annotations.get("solexa_quality")
        if q_old is not None and q_new is not None:
            # Mapping from Solexa to PHRED is lossy, but so is PHRED to Solexa.
            # Assume "old" is the original, and "new" has been converted.
            converted = [round(QualityIO.solexa_quality_from_phred(q)) for q in q_old]
            if truncate is not None:
                converted = [min(q, truncate) for q in converted]
            err_msg = f"mismatch converting phred_quality {q_old} to solexa_quality"
            if msg is not None:
                err_msg = f"{msg}: {err_msg}"
            self.assertEqual(converted, q_new, msg=err_msg)

        q_old = old.letter_annotations.get("solexa_quality")
        q_new = new.letter_annotations.get("phred_quality")
        if q_old is not None and q_new is not None:
            # Mapping from Solexa to PHRED is lossy, but so is PHRED to Solexa.
            # Assume "old" is the original, and "new" has been converted.
            converted = [round(QualityIO.phred_quality_from_solexa(q)) for q in q_old]
            if truncate is not None:
                converted = [min(q, truncate) for q in converted]
            err_msg = f"mismatch converting solexa_quality {q_old} to phred_quality"
            if msg is not None:
                err_msg = f"{msg}: {err_msg}"
            self.assertEqual(converted, q_new, msg=err_msg)


class TestFastqQualityDecoding(unittest.TestCase):
    """Test decoding across the complete valid range of each FASTQ variant."""

    def test_valid_quality_ranges(self):
        tests = [
            ("fastq-sanger", "phred_quality", range(94), 33),
            ("fastq-illumina", "phred_quality", range(63), 64),
            ("fastq-solexa", "solexa_quality", range(-5, 63), 64),
        ]
        for fmt, key, expected_range, offset in tests:
            with self.subTest(fmt=fmt):
                expected = list(expected_range)
                qualities = "".join(chr(score + offset) for score in expected)
                fastq = f"@test\n{'N' * len(expected)}\n+\n{qualities}\n"
                record = SeqIO.read(StringIO(fastq), fmt)
                decoded = record.letter_annotations[key]
                self.assertIs(type(decoded), list)
                self.assertEqual(decoded, expected)

    def test_invalid_quality_boundaries(self):
        tests = [
            ("fastq-sanger", " !", "\x7f!"),
            ("fastq-illumina", "?@", "\x7f@"),
            ("fastq-solexa", ":;", "\x7f;"),
        ]
        for fmt, below, above in tests:
            for qualities in (below, above):
                with self.subTest(fmt=fmt, qualities=qualities):
                    fastq = f"@test\nNN\n+\n{qualities}\n"
                    with self.assertRaises(QualityIO.InvalidCharError) as cm:
                        SeqIO.read(StringIO(fastq), fmt)
                    self.assertEqual(cm.exception.full_string, qualities)
                    self.assertEqual(cm.exception.index, 0)


class TestFastqCompactQualities(unittest.TestCase):
    """Test the opt-in compact=True array backing of FASTQ qualities."""

    iterators = {
        "fastq": QualityIO.FastqPhredIterator,
        "fastq-sanger": QualityIO.FastqPhredIterator,
        "fastq-illumina": QualityIO.FastqIlluminaIterator,
        "fastq-solexa": QualityIO.FastqSolexaIterator,
    }
    files = [
        (support.DATA / "Quality" / "example.fastq", "fastq", "phred_quality"),
        (
            support.DATA / "Quality" / "longreads_original_sanger.fastq",
            "fastq",
            "phred_quality",
        ),
        (
            support.DATA / "Quality" / "sanger_full_range_original_sanger.fastq",
            "fastq-sanger",
            "phred_quality",
        ),
        (
            support.DATA / "Quality" / "illumina_full_range_original_illumina.fastq",
            "fastq-illumina",
            "phred_quality",
        ),
        (
            support.DATA / "Quality" / "solexa_full_range_original_solexa.fastq",
            "fastq-solexa",
            "solexa_quality",
        ),
    ]

    def parse_both(self, path, fmt):
        """Parse a file twice, returning list-backed and array-backed records."""
        iterator = self.iterators[fmt]
        return list(iterator(path)), list(iterator(path, compact=True))

    def test_default_unchanged(self):
        for path, fmt, key in self.files:
            with self.subTest(path=path):
                for record in SeqIO.parse(path, fmt):
                    self.assertIs(type(record.letter_annotations[key]), list)
                for record in self.iterators[fmt](path):
                    self.assertIs(type(record.letter_annotations[key]), list)

    def test_valid_quality_ranges(self):
        tests = [
            ("fastq-sanger", "phred_quality", range(94), 33),
            ("fastq-illumina", "phred_quality", range(63), 64),
            ("fastq-solexa", "solexa_quality", range(-5, 63), 64),
        ]
        for fmt, key, expected_range, offset in tests:
            with self.subTest(fmt=fmt):
                expected = list(expected_range)
                qualities = "".join(chr(score + offset) for score in expected)
                fastq = f"@test\n{'N' * len(expected)}\n+\n{qualities}\n"
                iterator = self.iterators[fmt](StringIO(fastq), compact=True)
                decoded = next(iterator).letter_annotations[key]
                self.assertIs(type(decoded), array.array)
                self.assertEqual(decoded.typecode, "b")
                self.assertEqual(decoded.tolist(), expected)

    def test_invalid_still_rejected(self):
        for fmt in ("fastq-sanger", "fastq-illumina", "fastq-solexa"):
            with self.subTest(fmt=fmt):
                iterator = self.iterators[fmt](
                    StringIO("@test\nNN\n+\n\x7f!\n"), compact=True
                )
                with self.assertRaises(QualityIO.InvalidCharError):
                    next(iterator)

    def test_same_values_as_list(self):
        for path, fmt, key in self.files:
            with self.subTest(path=path):
                lists, arrays = self.parse_both(path, fmt)
                self.assertEqual(len(lists), len(arrays))
                for old, new in zip(lists, arrays):
                    self.assertEqual(old.id, new.id)
                    self.assertEqual(old.seq, new.seq)
                    self.assertEqual(
                        old.letter_annotations[key], list(new.letter_annotations[key])
                    )

    def test_write(self):
        formats = [
            "fastq",
            "fastq-sanger",
            "fastq-illumina",
            "fastq-solexa",
            "qual",
            "fasta",
        ]
        for path, fmt, key in self.files:
            lists, arrays = self.parse_both(path, fmt)
            for out_fmt in formats:
                with self.subTest(path=path, out_fmt=out_fmt):
                    outputs = []
                    for records in (lists, arrays):
                        handle = StringIO()
                        # Writing the full ranges to a narrower variant warns
                        # about truncation; both backings should warn alike.
                        with warnings.catch_warnings(record=True) as caught:
                            warnings.simplefilter("always")
                            SeqIO.write(records, handle, out_fmt)
                        messages = [str(w.message) for w in caught]
                        outputs.append((handle.getvalue(), messages))
                    self.assertEqual(outputs[0], outputs[1])

    def test_slicing(self):
        for path, fmt, key in self.files:
            with self.subTest(path=path):
                lists, arrays = self.parse_both(path, fmt)
                old, new = lists[0], arrays[0]
                for index in (
                    slice(3, 10),
                    slice(None, None, 2),
                    slice(None, None, -1),
                ):
                    sliced = new[index].letter_annotations[key]
                    self.assertIs(type(sliced), array.array)
                    self.assertEqual(list(sliced), old[index].letter_annotations[key])

    def test_reverse_complement(self):
        for path, fmt, key in self.files:
            with self.subTest(path=path):
                lists, arrays = self.parse_both(path, fmt)
                old, new = lists[0], arrays[0]
                rc = new.reverse_complement()
                self.assertIs(type(rc.letter_annotations[key]), array.array)
                self.assertEqual(
                    list(rc.letter_annotations[key]),
                    old.reverse_complement().letter_annotations[key],
                )

    def test_add(self):
        lists, arrays = self.parse_both(
            support.DATA / "Quality" / "example.fastq", "fastq"
        )
        combined = arrays[0] + arrays[1]
        qualities = combined.letter_annotations["phred_quality"]
        self.assertIs(type(qualities), array.array)
        self.assertEqual(
            list(qualities),
            lists[0].letter_annotations["phred_quality"]
            + lists[1].letter_annotations["phred_quality"],
        )
        # Mixing the two backings cannot concatenate, and says why:
        with self.assertRaisesRegex(TypeError, r"\(array and list\)"):
            arrays[0] + lists[1]

    def test_upper_lower(self):
        record = next(
            QualityIO.FastqPhredIterator(
                support.DATA / "Quality" / "example.fastq", compact=True
            )
        )
        before = list(record.letter_annotations["phred_quality"])
        for method in ("upper", "lower"):
            with self.subTest(method=method):
                derived = getattr(record, method)()
                qualities = derived.letter_annotations["phred_quality"]
                self.assertIs(type(qualities), array.array)
                self.assertEqual(list(qualities), before)
                # An independent copy, as for a list:
                qualities[0] = 0
                self.assertEqual(
                    list(record.letter_annotations["phred_quality"]), before
                )

    def test_letter_annotations_validation(self):
        record = next(
            QualityIO.FastqPhredIterator(
                support.DATA / "Quality" / "example.fastq", compact=True
            )
        )
        with self.assertRaises(TypeError):
            record.letter_annotations["phred_quality"] = array.array("b", [1, 2])
        record.letter_annotations["phred_quality"] = array.array(
            "b", range(len(record))
        )
        self.assertEqual(
            list(record.letter_annotations["phred_quality"]), list(range(len(record)))
        )

    def test_pickle(self):
        record = next(
            QualityIO.FastqSolexaIterator(
                support.DATA / "Quality" / "solexa_faked.fastq", compact=True
            )
        )
        clone = pickle.loads(pickle.dumps(record))
        self.assertEqual(clone.letter_annotations, record.letter_annotations)
        self.assertIs(type(clone.letter_annotations["solexa_quality"]), array.array)

    def test_index_get_raw(self):
        for path, fmt, key in self.files:
            index = SeqIO.index(path, fmt)
            self.addCleanup(index.close)
            with self.subTest(path=path):
                for name in index:
                    # SeqIO.index takes no parser options, so it stays list-backed,
                    # but get_raw hands back text the compact parser reads:
                    expected = index[name].letter_annotations[key]
                    self.assertIs(type(expected), list)
                    raw = index.get_raw(name).decode()
                    record = next(self.iterators[fmt](StringIO(raw), compact=True))
                    self.assertEqual(record.id, name)
                    self.assertEqual(list(record.letter_annotations[key]), expected)


class TestFastqErrors(unittest.TestCase):
    """Test reject invalid FASTQ files."""

    def check_fails(self, filename, good_count, formats=None, raw=True):
        if not formats:
            formats = ["fastq-sanger", "fastq-solexa", "fastq-illumina"]
        msg = f"SeqIO.parse failed to detect error in {filename}"
        for fmt in formats:
            records = SeqIO.parse(filename, fmt)
            for i in range(good_count):
                record = next(records)  # Make sure no errors!
                self.assertIsInstance(record, SeqRecord)
            # Detect error in the next record:
            with self.assertRaises(ValueError, msg=msg) as cm:
                record = next(records)

    def check_general_fails(self, filename, good_count):
        tuples = QualityIO.FastqGeneralIterator(filename)
        msg = f"FastqGeneralIterator failed to detect error in {filename}"
        for i in range(good_count):
            title, seq, qual = next(tuples)  # Make sure no errors!
        # Detect error in the next record:
        with self.assertRaises(ValueError, msg=msg) as cm:
            title, seq, qual = next(tuples)

    def check_general_passes(self, filename, record_count):
        tuples = QualityIO.FastqGeneralIterator(filename)
        # This "raw" parser doesn't check the ASCII characters which means
        # certain invalid FASTQ files will get parsed without errors.
        msg = f"FastqGeneralIterator failed to parse {filename}"
        count = 0
        for title, seq, qual in tuples:
            self.assertEqual(len(seq), len(qual), msg=msg)
            count += 1
        self.assertEqual(count, record_count, msg=msg)

    def test_reject_high_and_low(self):
        # These FASTQ files will be rejected by both the low level parser AND
        # the high level SeqRecord parser:
        tests = [
            (support.DATA / "Quality" / "error_diff_ids.fastq", 2),
            (support.DATA / "Quality" / "error_no_qual.fastq", 0),
            (support.DATA / "Quality" / "error_long_qual.fastq", 3),
            (support.DATA / "Quality" / "error_short_qual.fastq", 2),
            (support.DATA / "Quality" / "error_double_seq.fastq", 3),
            (support.DATA / "Quality" / "error_double_qual.fastq", 2),
            (support.DATA / "Quality" / "error_tabs.fastq", 0),
            (support.DATA / "Quality" / "error_spaces.fastq", 0),
            (support.DATA / "Quality" / "error_trunc_in_title.fastq", 4),
            (support.DATA / "Quality" / "error_trunc_in_seq.fastq", 4),
            (support.DATA / "Quality" / "error_trunc_in_plus.fastq", 4),
            (support.DATA / "Quality" / "error_trunc_in_qual.fastq", 4),
            (support.DATA / "Quality" / "error_trunc_at_seq.fastq", 4),
            (support.DATA / "Quality" / "error_trunc_at_plus.fastq", 4),
            (support.DATA / "Quality" / "error_trunc_at_qual.fastq", 4),
        ]
        for path, count in tests:
            self.check_fails(path, count)
            self.check_general_fails(path, count)

    def test_reject_high_but_not_low(self):
        # These FASTQ files which will be rejected by the high level SeqRecord
        # parser, but will be accepted by the low level parser:
        tests = [
            (support.DATA / "Quality" / "error_qual_del.fastq", 3, 5),
            (support.DATA / "Quality" / "error_qual_space.fastq", 3, 5),
            (support.DATA / "Quality" / "error_qual_vtab.fastq", 0, 5),
            (support.DATA / "Quality" / "error_qual_escape.fastq", 4, 5),
            (support.DATA / "Quality" / "error_qual_unit_sep.fastq", 2, 5),
            (support.DATA / "Quality" / "error_qual_tab.fastq", 4, 5),
            (support.DATA / "Quality" / "error_qual_null.fastq", 0, 5),
        ]
        for path, good_count, full_count in tests:
            self.check_fails(path, good_count)
            self.check_general_passes(path, full_count)


class TestReferenceSffConversions(unittest.TestCase):
    def check(self, sff_name, sff_format, out_name, fmt):
        wanted = list(SeqIO.parse(out_name, fmt))
        data = StringIO()
        count = SeqIO.convert(sff_name, sff_format, data, fmt)
        self.assertEqual(count, len(wanted))
        data.seek(0)
        converted = list(SeqIO.parse(data, fmt))
        self.assertEqual(len(wanted), len(converted))
        for old, new in zip(wanted, converted):
            self.assertEqual(old.id, new.id)
            self.assertEqual(old.name, new.name)
            if fmt != "qual":
                self.assertEqual(old.seq, new.seq)
            elif fmt != "fasta":
                self.assertEqual(
                    old.letter_annotations["phred_quality"],
                    new.letter_annotations["phred_quality"],
                )

    def check_sff(self, sff_name):
        self.check(
            sff_name,
            "sff",
            support.DATA / "Roche" / "E3MFGYR02_random_10_reads_no_trim.fasta",
            "fasta",
        )
        self.check(
            sff_name,
            "sff",
            support.DATA / "Roche" / "E3MFGYR02_random_10_reads_no_trim.qual",
            "qual",
        )
        self.check(
            sff_name,
            "sff-trim",
            support.DATA / "Roche" / "E3MFGYR02_random_10_reads.fasta",
            "fasta",
        )
        self.check(
            sff_name,
            "sff-trim",
            support.DATA / "Roche" / "E3MFGYR02_random_10_reads.qual",
            "qual",
        )

    def test_original(self):
        """Test converting E3MFGYR02_random_10_reads.sff into FASTA+QUAL."""
        self.check_sff(support.DATA / "Roche" / "E3MFGYR02_random_10_reads.sff")

    def test_no_manifest(self):
        """Test converting E3MFGYR02_no_manifest.sff into FASTA+QUAL."""
        self.check_sff(support.DATA / "Roche" / "E3MFGYR02_no_manifest.sff")

    def test_alt_index_at_start(self):
        """Test converting E3MFGYR02_alt_index_at_start into FASTA+QUAL."""
        self.check_sff(support.DATA / "Roche" / "E3MFGYR02_alt_index_at_start.sff")

    def test_alt_index_in_middle(self):
        """Test converting E3MFGYR02_alt_index_in_middle into FASTA+QUAL."""
        self.check_sff(support.DATA / "Roche" / "E3MFGYR02_alt_index_in_middle.sff")

    def test_alt_index_at_end(self):
        """Test converting E3MFGYR02_alt_index_at_end into FASTA+QUAL."""
        self.check_sff(support.DATA / "Roche" / "E3MFGYR02_alt_index_at_end.sff")

    def test_index_at_start(self):
        """Test converting E3MFGYR02_index_at_start into FASTA+QUAL."""
        self.check_sff(support.DATA / "Roche" / "E3MFGYR02_index_at_start.sff")

    def test_index_at_end(self):
        """Test converting E3MFGYR02_index_in_middle into FASTA+QUAL."""
        self.check_sff(support.DATA / "Roche" / "E3MFGYR02_index_in_middle.sff")


class TestReferenceFastqConversions(unittest.TestCase):
    """Tests where we have reference output."""

    def simple_check(self, base_name, in_variant):
        for out_variant in ["sanger", "solexa", "illumina"]:
            in_filename = (
                support.DATA / "Quality" / f"{base_name}_original_{in_variant}.fastq"
            )
            self.assertTrue(os.path.isfile(in_filename))
            # Load the reference output...
            with open(
                support.DATA / "Quality" / f"{base_name}_as_{out_variant}.fastq"
            ) as handle:
                expected = handle.read()

            with warnings.catch_warnings():
                if out_variant != "sanger":
                    # Ignore data loss warnings from max qualities
                    warnings.simplefilter("ignore", BiopythonWarning)
                # Check matches using convert...
                handle = StringIO()
                SeqIO.convert(
                    in_filename, "fastq-" + in_variant, handle, "fastq-" + out_variant
                )
                self.assertEqual(expected, handle.getvalue())
                # Check matches using parse/write
                handle = StringIO()
                SeqIO.write(
                    SeqIO.parse(in_filename, "fastq-" + in_variant),
                    handle,
                    "fastq-" + out_variant,
                )
                self.assertEqual(expected, handle.getvalue())

    def test_reference_conversion(self):
        tests = [
            ("illumina_full_range", "illumina"),
            ("sanger_full_range", "sanger"),
            ("longreads", "sanger"),
            ("solexa_full_range", "solexa"),
            ("misc_dna", "sanger"),
            ("wrapping", "sanger"),
            ("misc_rna", "sanger"),
        ]
        for base_name, variant in tests:
            assert variant in ["sanger", "solexa", "illumina"]
            self.simple_check(base_name, variant)


class TestQual(QualityIOTestBaseClass):
    """Tests with QUAL files."""

    def test_paired(self):
        """Check FASTQ parsing matches FASTA+QUAL parsing."""
        with (
            open(support.DATA / "Quality" / "example.fasta") as f,
            open(support.DATA / "Quality" / "example.qual") as q,
        ):
            records1 = list(QualityIO.PairedFastaQualIterator(f, q))
        records2 = list(
            SeqIO.parse(support.DATA / "Quality" / "example.fastq", "fastq")
        )
        self.compare_records(records1, records2)

    def test_qual_empty_title(self):
        """Check QUAL record with empty title."""
        record = SeqIO.read(StringIO(">\n1 2 3"), "qual")
        self.assertEqual(record.id, "")
        self.assertEqual(record.description, "")
        self.assertEqual(record.letter_annotations["phred_quality"], [1, 2, 3])

    def test_qual(self):
        """Check FASTQ parsing matches QUAL parsing."""
        records1 = list(SeqIO.parse(support.DATA / "Quality" / "example.qual", "qual"))
        records2 = list(
            SeqIO.parse(support.DATA / "Quality" / "example.fastq", "fastq")
        )
        # Will ignore the unknown sequences :)
        self.compare_records(records1, records2)

    def test_qual_out(self):
        """Check FASTQ to QUAL output."""
        records = SeqIO.parse(support.DATA / "Quality" / "example.fastq", "fastq")
        h = StringIO()
        SeqIO.write(records, h, "qual")
        with open(support.DATA / "Quality" / "example.qual") as expected:
            self.assertEqual(h.getvalue(), expected.read())

    def test_fasta(self):
        """Check FASTQ parsing matches FASTA parsing."""
        records1 = list(
            SeqIO.parse(support.DATA / "Quality" / "example.fasta", "fasta")
        )
        records2 = list(
            SeqIO.parse(support.DATA / "Quality" / "example.fastq", "fastq")
        )
        self.compare_records(records1, records2)

    def test_fasta_out(self):
        """Check FASTQ to FASTA output."""
        records = SeqIO.parse(support.DATA / "Quality" / "example.fastq", "fastq")
        h = StringIO()
        SeqIO.write(records, h, "fasta")
        with open(support.DATA / "Quality" / "example.fasta") as expected:
            self.assertEqual(h.getvalue(), expected.read())

    def test_qual_negative(self):
        """Check QUAL negative scores mapped to PHRED zero."""
        data = """>1117_10_107_F3
23 31 -1 -1 -1 29 -1 -1 20 32 -1 18 25 7 -1 6 -1 -1 -1 30 -1 20 13 7 -1 -1 21 30 -1 24 -1 22 -1 -1 22 14 -1 12 26 21 -1 5 -1 -1 -1 20 -1 -1 12 28
>1117_10_146_F3
20 33 -1 -1 -1 29 -1 -1 28 28 -1 7 16 5 -1 30 -1 -1 -1 14 -1 4 13 4 -1 -1 11 13 -1 5 -1 7 -1 -1 10 16 -1 4 12 15 -1 8 -1 -1 -1 16 -1 -1 10 4
>1117_10_1017_F3
33 33 -1 -1 -1 27 -1 -1 17 16 -1 28 24 11 -1 6 -1 -1 -1 29 -1 8 29 24 -1 -1 8 8 -1 20 -1 13 -1 -1 8 13 -1 28 10 24 -1 10 -1 -1 -1 4 -1 -1 7 6
>1117_11_136_F3
16 22 -1 -1 -1 33 -1 -1 30 27 -1 27 28 32 -1 29 -1 -1 -1 27 -1 18 9 6 -1 -1 23 16 -1 26 -1 5 7 -1 22 7 -1 18 14 8 -1 8 -1 -1 -1 11 -1 -1 4 24"""  # noqa : W291
        h = StringIO(data)
        h2 = StringIO()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonParserWarning)
            records = SeqIO.parse(h, "qual")

            def add_sequence(records):
                for record in records:
                    record.seq = Seq(len(record.seq) * "?")
                    yield record

            records = add_sequence(records)
            self.assertEqual(4, SeqIO.write(records, h2, "fastq"))
        self.assertEqual(
            h2.getvalue(),
            """\
@1117_10_107_F3
??????????????????????????????????????????????????
+
8@!!!>!!5A!3:(!'!!!?!5.(!!6?!9!7!!7/!-;6!&!!!5!!-=
@1117_10_146_F3
??????????????????????????????????????????????????
+
5B!!!>!!==!(1&!?!!!/!%.%!!,.!&!(!!+1!%-0!)!!!1!!+%
@1117_10_1017_F3
??????????????????????????????????????????????????
+
BB!!!<!!21!=9,!'!!!>!)>9!!))!5!.!!).!=+9!+!!!%!!('
@1117_11_136_F3
??????????????????????????????????????????????????
+
17!!!B!!?<!<=A!>!!!<!3*'!!81!;!&(!7(!3/)!)!!!,!!%9
""",
        )


class TestReadWrite(unittest.TestCase):
    """Test can read and write back files."""

    def test_fastq_2000(self):
        """Read and write back simple example with upper case 2000bp read."""
        data = f"@{'id descr goes here'}\n{'ACGT' * 500}\n+\n{'!@a~' * 500}\n"
        handle = StringIO()
        self.assertEqual(
            1, SeqIO.write(SeqIO.parse(StringIO(data), "fastq"), handle, "fastq")
        )
        self.assertEqual(data, handle.getvalue())

    def test_fastq_1000(self):
        """Read and write back simple example with mixed case 1000bp read."""
        data = "@%s\n%s\n+\n%s\n" % (
            "id descr goes here",
            "ACGTNncgta" * 100,
            "abcd!!efgh" * 100,
        )
        handle = StringIO()
        self.assertEqual(
            1, SeqIO.write(SeqIO.parse(StringIO(data), "fastq"), handle, "fastq")
        )
        self.assertEqual(data, handle.getvalue())

    def test_fastq_dna(self):
        """Read and write back simple example with ambiguous DNA."""
        # First in upper case...
        data = "@%s\n%s\n+\n%s\n" % (
            "id descr goes here",
            ambiguous_dna_letters.upper(),
            "".join(chr(33 + q) for q in range(len(ambiguous_dna_letters))),
        )
        handle = StringIO()
        self.assertEqual(
            1, SeqIO.write(SeqIO.parse(StringIO(data), "fastq"), handle, "fastq")
        )
        self.assertEqual(data, handle.getvalue())
        # Now in lower case...
        data = "@%s\n%s\n+\n%s\n" % (
            "id descr goes here",
            ambiguous_dna_letters.lower(),
            "".join(chr(33 + q) for q in range(len(ambiguous_dna_letters))),
        )
        handle = StringIO()
        self.assertEqual(
            1, SeqIO.write(SeqIO.parse(StringIO(data), "fastq"), handle, "fastq")
        )
        self.assertEqual(data, handle.getvalue())

    def test_fastq_rna(self):
        """Read and write back simple example with ambiguous RNA."""
        # First in upper case...
        data = "@%s\n%s\n+\n%s\n" % (
            "id descr goes here",
            ambiguous_rna_letters.upper(),
            "".join(chr(33 + q) for q in range(len(ambiguous_rna_letters))),
        )
        handle = StringIO()
        self.assertEqual(
            1, SeqIO.write(SeqIO.parse(StringIO(data), "fastq"), handle, "fastq")
        )
        self.assertEqual(data, handle.getvalue())
        # Now in lower case...
        data = "@%s\n%s\n+\n%s\n" % (
            "id descr goes here",
            ambiguous_rna_letters.lower(),
            "".join(chr(33 + q) for q in range(len(ambiguous_rna_letters))),
        )
        handle = StringIO()
        self.assertEqual(
            1, SeqIO.write(SeqIO.parse(StringIO(data), "fastq"), handle, "fastq")
        )
        self.assertEqual(data, handle.getvalue())


class TestWriteRead(QualityIOTestBaseClass):
    """Test can write and read back files."""

    def write_read(self, filename, in_format, out_format):
        records = list(SeqIO.parse(filename, in_format))
        for record in records:
            try:
                bytes(record.seq)
            except UndefinedSequenceError:
                record.seq = Seq("N" * len(record.seq))
        mode = self.get_mode(out_format)
        if mode == "b":
            handle = BytesIO()
        else:
            handle = StringIO()
        SeqIO.write(records, handle, out_format)
        handle.seek(0)
        # Now load it back and check it agrees,
        records2 = list(SeqIO.parse(handle, out_format))
        self.compare_records(records, records2, out_format)

    def test_generated(self):
        """Write and read back odd SeqRecord objects."""
        record1 = SeqRecord(
            Seq("ACGT" * 500),
            id="Test",
            description="Long " * 500,
            letter_annotations={"phred_quality": [40, 30, 20, 10] * 500},
        )
        record2 = SeqRecord(
            MutableSeq("NGGC" * 1000),
            id="Mut",
            description="very " * 1000 + "long",
            letter_annotations={"phred_quality": [0, 5, 5, 10] * 1000},
        )
        record3 = SeqRecord(
            Seq("N" * 2000),
            id="Unk",
            description="l" + ("o" * 1000) + "ng",
            letter_annotations={"phred_quality": [0, 1] * 1000},
        )
        record4 = SeqRecord(
            Seq("ACGT" * 500),
            id="no_descr",
            description="",
            name="",
            letter_annotations={"phred_quality": [40, 50, 60, 62] * 500},
        )
        record5 = SeqRecord(
            Seq(""),
            id="empty_p",
            description="(could have been trimmed lots)",
            letter_annotations={"phred_quality": []},
        )
        record6 = SeqRecord(
            Seq(""),
            id="empty_s",
            description="(could have been trimmed lots)",
            letter_annotations={"solexa_quality": []},
        )
        record7 = SeqRecord(
            Seq("ACNN" * 500),
            id="Test_Sol",
            description="Long " * 500,
            letter_annotations={"solexa_quality": [40, 30, 0, -5] * 500},
        )
        record8 = SeqRecord(
            Seq("ACGT"),
            id="HighQual",
            description="With very large qualities that even Sanger FASTQ can't hold!",
            letter_annotations={"solexa_quality": [0, 10, 100, 1000]},
        )
        # TODO - Record with no identifier?
        records = [
            record1,
            record2,
            record3,
            record4,
            record5,
            record6,
            record7,
            record8,
        ]
        for fmt in ["fasta", "fastq", "fastq-solexa", "fastq-illumina", "qual"]:
            handle = StringIO()
            with warnings.catch_warnings():
                # TODO - Have a Biopython defined "DataLossWarning?"
                warnings.simplefilter("ignore", BiopythonWarning)
                SeqIO.write(records, handle, fmt)
            handle.seek(0)
            self.compare_records(records, list(SeqIO.parse(handle, fmt)), fmt)

    def check(self, filename, fmt, out_formats):
        for f in out_formats:
            self.write_read(filename, fmt, f)

    def test_tricky(self):
        """Write and read back tricky.fastq."""
        self.check(
            support.DATA / "Quality" / "tricky.fastq",
            "fastq",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
            ],
        )

    def test_sanger_93(self):
        """Write and read back sanger_93.fastq."""
        self.check(
            support.DATA / "Quality" / "sanger_93.fastq",
            "fastq",
            ["fastq", "fastq-sanger", "fasta", "qual", "phd"],
        )
        with warnings.catch_warnings():
            # TODO - Have a Biopython defined "DataLossWarning?"
            warnings.simplefilter("ignore", BiopythonWarning)
            self.check(
                support.DATA / "Quality" / "sanger_93.fastq",
                "fastq",
                ["fastq-solexa", "fastq-illumina"],
            )

    def test_sanger_faked(self):
        """Write and read back sanger_faked.fastq."""
        self.check(
            support.DATA / "Quality" / "sanger_faked.fastq",
            "fastq",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
            ],
        )

    def test_example_fasta(self):
        """Write and read back example.fasta."""
        self.write_read(support.DATA / "Quality" / "example.fasta", "fasta", "fasta")
        # TODO - tests to check can't write FASTQ or QUAL...

    def test_example_fastq(self):
        """Write and read back example.fastq."""
        self.check(
            support.DATA / "Quality" / "example.fastq",
            "fastq",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
            ],
        )

    def test_example_qual(self):
        """Write and read back example.qual."""
        self.check(
            support.DATA / "Quality" / "example.qual",
            "qual",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
            ],
        )

    def test_solexa_faked(self):
        """Write and read back solexa_faked.fastq."""
        self.check(
            support.DATA / "Quality" / "solexa_faked.fastq",
            "fastq-solexa",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
            ],
        )

    def test_solexa_example(self):
        """Write and read back solexa_example.fastq."""
        self.check(
            support.DATA / "Quality" / "solexa_example.fastq",
            "fastq-solexa",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
            ],
        )

    def test_illumina_faked(self):
        """Write and read back illumina_faked.fastq."""
        self.check(
            support.DATA / "Quality" / "illumina_faked.fastq",
            "fastq-illumina",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
            ],
        )

    def test_greek_sff(self):
        """Write and read back greek.sff."""
        self.check(
            support.DATA / "Roche" / "greek.sff",
            "sff",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
                "sff",
            ],
        )

    def test_paired_sff(self):
        """Write and read back paired.sff."""
        self.check(
            support.DATA / "Roche" / "paired.sff",
            "sff",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
                "sff",
            ],
        )

    def test_E3MFGYR02(self):
        """Write and read back E3MFGYR02_random_10_reads.sff."""
        self.check(
            support.DATA / "Roche" / "E3MFGYR02_random_10_reads.sff",
            "sff",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
                "sff",
            ],
        )

    def test_E3MFGYR02_no_manifest(self):
        """Write and read back E3MFGYR02_no_manifest.sff."""
        self.check(
            support.DATA / "Roche" / "E3MFGYR02_no_manifest.sff",
            "sff",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
                "sff",
            ],
        )

    def test_E3MFGYR02_index_at_start(self):
        """Write and read back E3MFGYR02_index_at_start.sff."""
        self.check(
            support.DATA / "Roche" / "E3MFGYR02_index_at_start.sff",
            "sff",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
                "sff",
            ],
        )

    def test_E3MFGYR02_index_in_middle(self):
        """Write and read back E3MFGYR02_index_in_middle.sff."""
        self.check(
            support.DATA / "Roche" / "E3MFGYR02_index_in_middle.sff",
            "sff",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
                "sff",
            ],
        )

    def test_E3MFGYR02_alt_index_at_start(self):
        """Write and read back E3MFGYR02_alt_index_at_start.sff."""
        self.check(
            support.DATA / "Roche" / "E3MFGYR02_alt_index_at_start.sff",
            "sff",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
                "sff",
            ],
        )

    def test_E3MFGYR02_alt_index_in_middle(self):
        """Write and read back E3MFGYR02_alt_index_in_middle.sff."""
        self.check(
            support.DATA / "Roche" / "E3MFGYR02_alt_index_in_middle.sff",
            "sff",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
                "sff",
            ],
        )

    def test_E3MFGYR02_alt_index_at_end(self):
        """Write and read back E3MFGYR02_alt_index_at_end.sff."""
        self.check(
            support.DATA / "Roche" / "E3MFGYR02_alt_index_at_end.sff",
            "sff",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
                "sff",
            ],
        )

    def test_E3MFGYR02_trimmed(self):
        """Write and read back E3MFGYR02_random_10_reads.sff (trimmed)."""
        self.check(
            support.DATA / "Roche" / "E3MFGYR02_random_10_reads.sff",
            "sff-trim",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
            ],
        )  # not sff as output


class MappingTests(unittest.TestCase):
    """Quality mapping tests."""

    def test_solexa_quality_from_phred(self):
        """Mapping check for function solexa_quality_from_phred."""
        self.assertEqual(-5, round(QualityIO.solexa_quality_from_phred(0)))
        self.assertEqual(-5, round(QualityIO.solexa_quality_from_phred(1)))
        self.assertEqual(-2, round(QualityIO.solexa_quality_from_phred(2)))
        self.assertEqual(0, round(QualityIO.solexa_quality_from_phred(3)))
        self.assertEqual(2, round(QualityIO.solexa_quality_from_phred(4)))
        self.assertEqual(3, round(QualityIO.solexa_quality_from_phred(5)))
        self.assertEqual(5, round(QualityIO.solexa_quality_from_phred(6)))
        self.assertEqual(6, round(QualityIO.solexa_quality_from_phred(7)))
        self.assertEqual(7, round(QualityIO.solexa_quality_from_phred(8)))
        self.assertEqual(8, round(QualityIO.solexa_quality_from_phred(9)))
        for i in range(10, 100):
            self.assertEqual(i, round(QualityIO.solexa_quality_from_phred(i)))

    def test_phred_quality_from_solexa(self):
        """Mapping check for function phred_quality_from_solexa."""
        self.assertEqual(1, round(QualityIO.phred_quality_from_solexa(-5)))
        self.assertEqual(1, round(QualityIO.phred_quality_from_solexa(-4)))
        self.assertEqual(2, round(QualityIO.phred_quality_from_solexa(-3)))
        self.assertEqual(2, round(QualityIO.phred_quality_from_solexa(-2)))
        self.assertEqual(3, round(QualityIO.phred_quality_from_solexa(-1)))
        self.assertEqual(3, round(QualityIO.phred_quality_from_solexa(0)))
        self.assertEqual(4, round(QualityIO.phred_quality_from_solexa(1)))
        self.assertEqual(4, round(QualityIO.phred_quality_from_solexa(2)))
        self.assertEqual(5, round(QualityIO.phred_quality_from_solexa(3)))
        self.assertEqual(5, round(QualityIO.phred_quality_from_solexa(4)))
        self.assertEqual(6, round(QualityIO.phred_quality_from_solexa(5)))
        self.assertEqual(7, round(QualityIO.phred_quality_from_solexa(6)))
        self.assertEqual(8, round(QualityIO.phred_quality_from_solexa(7)))
        self.assertEqual(9, round(QualityIO.phred_quality_from_solexa(8)))
        self.assertEqual(10, round(QualityIO.phred_quality_from_solexa(9)))
        for i in range(10, 100):
            self.assertEqual(i, round(QualityIO.phred_quality_from_solexa(i)))

    def test_sanger_to_solexa(self):
        """Mapping check for FASTQ Sanger (0 to 93) to Solexa (-5 to 62)."""
        # The point of this test is the writing code doesn't actually use the
        # solexa_quality_from_phred function directly. For speed it uses a
        # cached dictionary of the mappings.
        seq = "N" * 94
        qual = "".join(chr(33 + q) for q in range(94))
        expected_sol = [
            min(62, int(round(QualityIO.solexa_quality_from_phred(q))))
            for q in range(94)
        ]
        in_handle = StringIO(f"@Test\n{seq}\n+\n{qual}")
        out_handle = StringIO()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", BiopythonWarning)
            SeqIO.write(
                SeqIO.parse(in_handle, "fastq-sanger"), out_handle, "fastq-solexa"
            )
            self.assertLessEqual(len(w), 1, w)
        out_handle.seek(0)
        record = SeqIO.read(out_handle, "fastq-solexa")
        self.assertEqual(record.seq, seq)
        self.assertEqual(record.letter_annotations["solexa_quality"], expected_sol)

    def test_solexa_to_sanger(self):
        """Mapping check for FASTQ Solexa (-5 to 62) to Sanger (0 to 62)."""
        # The point of this test is the writing code doesn't actually use the
        # solexa_quality_from_phred function directly. For speed it uses a
        # cached dictionary of the mappings.
        seq = "N" * 68
        qual = "".join(chr(64 + q) for q in range(-5, 63))
        expected_phred = [
            round(QualityIO.phred_quality_from_solexa(q)) for q in range(-5, 63)
        ]
        in_handle = StringIO(f"@Test\n{seq}\n+\n{qual}")
        out_handle = StringIO()
        SeqIO.write(SeqIO.parse(in_handle, "fastq-solexa"), out_handle, "fastq-sanger")
        out_handle.seek(0)
        record = SeqIO.read(out_handle, "fastq-sanger")
        self.assertEqual(record.seq, seq)
        self.assertEqual(record.letter_annotations["phred_quality"], expected_phred)

    def test_sanger_to_illumina(self):
        """Mapping check for FASTQ Sanger (0 to 93) to Illumina (0 to 62)."""
        seq = "N" * 94
        qual = "".join(chr(33 + q) for q in range(94))
        expected_phred = [min(62, q) for q in range(94)]
        in_handle = StringIO(f"@Test\n{seq}\n+\n{qual}")
        out_handle = StringIO()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", BiopythonWarning)
            SeqIO.write(
                SeqIO.parse(in_handle, "fastq-sanger"), out_handle, "fastq-illumina"
            )
            self.assertLessEqual(len(w), 1, w)
        out_handle.seek(0)
        record = SeqIO.read(out_handle, "fastq-illumina")
        self.assertEqual(record.seq, seq)
        self.assertEqual(record.letter_annotations["phred_quality"], expected_phred)

    def test_illumina_to_sanger(self):
        """Mapping check for FASTQ Illumina (0 to 62) to Sanger (0 to 62)."""
        seq = "N" * 63
        qual = "".join(chr(64 + q) for q in range(63))
        expected_phred = list(range(63))
        in_handle = StringIO(f"@Test\n{seq}\n+\n{qual}")
        out_handle = StringIO()
        SeqIO.write(
            SeqIO.parse(in_handle, "fastq-illumina"), out_handle, "fastq-sanger"
        )
        out_handle.seek(0)
        record = SeqIO.read(out_handle, "fastq-sanger")
        self.assertEqual(record.seq, seq)
        self.assertEqual(record.letter_annotations["phred_quality"], expected_phred)


class TestSFF(unittest.TestCase):
    """Test SFF specific details."""

    def test_overlapping_clip(self):
        record = next(SeqIO.parse(support.DATA / "Roche" / "greek.sff", "sff"))
        self.assertEqual(len(record), 395)
        s = record.seq.lower()
        # Apply overlapping clipping
        record.annotations["clip_qual_left"] = 51
        record.annotations["clip_qual_right"] = 44
        record.annotations["clip_adapter_left"] = 50
        record.annotations["clip_adapter_right"] = 75
        self.assertEqual(len(record), 395)
        self.assertEqual(len(record.seq), 395)
        # Save the clipped record...
        h = BytesIO()
        count = SeqIO.write(record, h, "sff")
        self.assertEqual(count, 1)
        # Now reload it...
        h.seek(0)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", BiopythonParserWarning)
            record = SeqIO.read(h, "sff")
            self.assertEqual(len(w), 1, w)
        self.assertEqual(record.annotations["clip_qual_left"], 51)
        self.assertEqual(record.annotations["clip_qual_right"], 44)
        self.assertEqual(record.annotations["clip_adapter_left"], 50)
        self.assertEqual(record.annotations["clip_adapter_right"], 75)
        self.assertEqual(len(record), 395)
        self.assertEqual(s, record.seq.lower())
        # And check with trimming applied...
        h.seek(0)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", BiopythonParserWarning)
            record = SeqIO.read(h, "sff-trim")
            self.assertEqual(len(w), 1, w)
        self.assertEqual(len(record), 0)

    def test_negative_clip(self):
        for clip in [
            "clip_qual_left",
            "clip_qual_right",
            "clip_adapter_left",
            "clip_adapter_right",
        ]:
            record = next(SeqIO.parse(support.DATA / "Roche" / "greek.sff", "sff"))
            self.assertEqual(len(record), 395)
            self.assertLessEqual(0, record.annotations[clip])
            record.annotations[clip] = -1
            with BytesIO() as h:
                self.assertRaises(ValueError, SeqIO.write, record, h, "sff")


class NonFastqTests(unittest.TestCase):
    def test_fasta_as_fastq(self):
        for f in ("fastq", "fastq-sanger", "fastq-solexa", "fastq-illumina"):
            generator = SeqIO.parse(support.DATA / "Fasta" / "elderberry.nu", f)
            self.assertRaises(ValueError, next, generator)

    def test_sff_as_fastq(self):
        for f in ("fastq", "fastq-sanger", "fastq-solexa", "fastq-illumina"):
            generator = SeqIO.parse(support.DATA / "Roche" / "greek.sff", f)
            self.assertRaises(ValueError, next, generator)


class TestsConverter(SeqIOConverterTestBaseClass, QualityIOTestBaseClass):
    def check_conversion(self, filename, in_format, out_format):
        msg = f"Convert {filename} from {in_format} to {out_format}"
        records = list(SeqIO.parse(filename, in_format))
        # Write it out...
        handle = StringIO()
        with warnings.catch_warnings():
            if out_format in (
                "fastq",
                "fastq-sanger",
                "fastq-solexa",
                "fastq-illumina",
            ):
                warnings.simplefilter("ignore", BiopythonWarning)
            SeqIO.write(records, handle, out_format)
        handle.seek(0)
        # Now load it back and check it agrees,
        records2 = list(SeqIO.parse(handle, out_format))
        self.assertEqual(len(records), len(records2), msg=msg)
        for record1, record2 in zip(records, records2):
            self.compare_record(record1, record2, out_format, msg=msg)
        # Finally, use the convert function, and check that agrees:
        handle2 = StringIO()
        with warnings.catch_warnings():
            if out_format in (
                "fastq",
                "fastq-sanger",
                "fastq-solexa",
                "fastq-illumina",
            ):
                warnings.simplefilter("ignore", BiopythonWarning)
            SeqIO.convert(filename, in_format, handle2, out_format)
        # We could re-parse this, but it is simpler and stricter:
        self.assertEqual(handle.getvalue(), handle2.getvalue(), msg=msg)

    def failure_check(self, filename, in_format, out_format):
        msg = "Confirm failure detection converting %s from %s to %s" % (
            filename,
            in_format,
            out_format,
        )
        # We want the SAME error message from parse/write as convert!
        with self.assertRaises(ValueError, msg=msg) as cm:
            records = list(SeqIO.parse(filename, in_format))
            self.write_records(records, out_format)
        err1 = str(cm.exception)
        # Now do the conversion...
        with self.assertRaises(ValueError, msg=msg) as cm:
            handle = StringIO()
            SeqIO.convert(filename, in_format, handle, out_format)
        err2 = str(cm.exception)
        # Verify that parse/write and convert give the same failure
        err_msg = f"{msg}: parse/write and convert gave different failures"
        self.assertEqual(err1, err2, msg=err_msg)

    def test_conversion(self):
        tests = [
            (support.DATA / "Quality" / "example.fastq", "fastq"),
            (support.DATA / "Quality" / "example.fastq", "fastq-sanger"),
            (support.DATA / "Quality" / "tricky.fastq", "fastq"),
            (support.DATA / "Quality" / "sanger_93.fastq", "fastq-sanger"),
            (support.DATA / "Quality" / "sanger_faked.fastq", "fastq-sanger"),
            (support.DATA / "Quality" / "solexa_faked.fastq", "fastq-solexa"),
            (support.DATA / "Quality" / "illumina_faked.fastq", "fastq-illumina"),
        ]
        for filename, fmt in tests:
            for in_format, out_format in self.formats:
                if in_format != fmt:
                    continue
                self.check_conversion(filename, in_format, out_format)

    def test_failure_detection(self):
        tests = [
            (support.DATA / "Quality" / "error_diff_ids.fastq", "fastq"),
            (support.DATA / "Quality" / "error_long_qual.fastq", "fastq"),
            (support.DATA / "Quality" / "error_no_qual.fastq", "fastq"),
            (support.DATA / "Quality" / "error_qual_del.fastq", "fastq"),
            (support.DATA / "Quality" / "error_qual_escape.fastq", "fastq"),
            (support.DATA / "Quality" / "error_qual_null.fastq", "fastq"),
            (support.DATA / "Quality" / "error_qual_space.fastq", "fastq"),
            (support.DATA / "Quality" / "error_qual_tab.fastq", "fastq"),
            (support.DATA / "Quality" / "error_qual_unit_sep.fastq", "fastq"),
            (support.DATA / "Quality" / "error_qual_vtab.fastq", "fastq"),
            (support.DATA / "Quality" / "error_short_qual.fastq", "fastq"),
            (support.DATA / "Quality" / "error_spaces.fastq", "fastq"),
            (support.DATA / "Quality" / "error_tabs.fastq", "fastq"),
            (support.DATA / "Quality" / "error_trunc_at_plus.fastq", "fastq"),
            (support.DATA / "Quality" / "error_trunc_at_qual.fastq", "fastq"),
            (support.DATA / "Quality" / "error_trunc_at_seq.fastq", "fastq"),
            (support.DATA / "Quality" / "error_trunc_in_title.fastq", "fastq"),
            (support.DATA / "Quality" / "error_trunc_in_seq.fastq", "fastq"),
            (support.DATA / "Quality" / "error_trunc_in_plus.fastq", "fastq"),
            (support.DATA / "Quality" / "error_trunc_in_qual.fastq", "fastq"),
            (support.DATA / "Quality" / "error_double_seq.fastq", "fastq"),
            (support.DATA / "Quality" / "error_double_qual.fastq", "fastq"),
        ]
        for filename, fmt in tests:
            for in_format, out_format in self.formats:
                if in_format != fmt:
                    continue
                if (
                    in_format
                    in ["fastq", "fastq-sanger", "fastq-solexa", "fastq-illumina"]
                    and out_format in ["fasta", "tab"]
                    and filename.name.startswith("error_qual_")
                ):
                    # TODO? These conversions don't check for bad characters in the quality,
                    # and in order to pass this strict test they should.
                    continue
                self.failure_check(filename, in_format, out_format)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
