# Copyright 2012 by Jeff Hussmann.  All rights reserved.
# Revisions copyright 2013-2016 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SffIO module."""

import re
import unittest
import warnings
from io import BytesIO

import support

from Bio import BiopythonParserWarning
from Bio import SeqIO
from Bio.SeqIO.SffIO import _sff_do_slow_index
from Bio.SeqIO.SffIO import _sff_find_roche_index
from Bio.SeqIO.SffIO import _sff_read_raw_record
from Bio.SeqIO.SffIO import _sff_read_roche_index
from Bio.SeqIO.SffIO import ReadRocheXmlManifest
from Bio.SeqIO.SffIO import SffIterator
from Bio.SeqIO.SffIO import SffWriter

# sffinfo E3MFGYR02_random_10_reads.sff | sed -n '/>\|Run Prefix\|Region\|XY/p'
test_data = """
>E3MFGYR02JWQ7T
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  3946_2103
>E3MFGYR02JA6IL
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  3700_3115
>E3MFGYR02JHD4H
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  3771_2095
>E3MFGYR02GFKUC
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2520_2738
>E3MFGYR02FTGED
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2268_2739
>E3MFGYR02FR9G7
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2255_0361
>E3MFGYR02GAZMS
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2468_1618
>E3MFGYR02HHZ8O
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2958_1574
>E3MFGYR02GPGB1
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2633_0607
>E3MFGYR02F7Z7G
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2434_1658"""


class TestUAN(unittest.TestCase):
    """Test annotations."""

    def setUp(self):
        self.records = list(
            SeqIO.parse(support.DATA / "Roche" / "E3MFGYR02_random_10_reads.sff", "sff")
        )
        self.test_annotations = {}
        for line in test_data.splitlines():
            fields = re.split(r"\s+", line.strip())
            if ">" in line:
                current_name = fields[0].lstrip(">")
                self.test_annotations[current_name] = {}
            elif "Prefix" in line:
                time_list = [int(v) for v in fields[2].split("_")[1:-1]]
                self.test_annotations[current_name]["time"] = time_list
            elif "Region" in line:
                region = int(fields[-1])
                self.test_annotations[current_name]["region"] = region
            elif "XY" in line:
                x, y = (int(v) for v in fields[-1].split("_"))
                self.test_annotations[current_name]["coords"] = (x, y)

    def test_time(self):
        for record in self.records:
            self.assertEqual(
                record.annotations["time"], self.test_annotations[record.name]["time"]
            )

    def test_region(self):
        for record in self.records:
            self.assertEqual(
                record.annotations["region"],
                self.test_annotations[record.name]["region"],
            )

    def test_coords(self):
        for record in self.records:
            self.assertEqual(
                record.annotations["coords"],
                self.test_annotations[record.name]["coords"],
            )


class TestErrors(unittest.TestCase):
    with open(support.DATA / "Roche" / "E3MFGYR02_random_10_reads.sff", "rb") as handle:
        good = handle.read()

    def test_empty(self):
        fh = BytesIO()
        with self.assertRaises(ValueError) as cm:
            records = list(SeqIO.parse(fh, "sff"))
        self.assertEqual(str(cm.exception), "Empty file.")

    def check_bad_header(self, header, msg):
        with self.assertRaises(ValueError) as cm:
            records = list(SeqIO.parse(BytesIO(header), "sff"))
        err = str(cm.exception)
        if isinstance(msg, (tuple, list)):
            self.assertIn(err, msg, f"Unexpected error: {err}")
        else:
            self.assertEqual(err, msg)

    def test_30bytes(self):
        self.check_bad_header(b"x" * 30, "File too small to hold a valid SFF header.")

    def test_31bytes(self):
        self.check_bad_header(
            b"x" * 31,
            (
                "SFF file did not start '.sff', but 'xxxx'",
                "SFF file did not start '.sff', but b'xxxx'",
            ),
        )

    def test_31bytes_index_header(self):
        self.check_bad_header(
            b".srt" + b"x" * 27, "Handle seems to be at SFF index block, not start"
        )

    def test_31bytes_bad_ver(self):
        self.check_bad_header(
            b".sff1.00" + b"x" * 23, "Unsupported SFF version in header, 49.46.48.48"
        )

    def test_31bytes_bad_flowgram(self):
        self.check_bad_header(
            b".sff\x00\x00\x00\x01" + b"x" * 23,
            "Flowgram format code 120 not supported",
        )

    def test_bad_index_offset(self):
        bad = self.good[:12] + b"\x00\x00\x00\x00" + self.good[16:]
        self.check_bad_header(bad, "Index offset 0 but index length 764")

    def test_bad_index_length(self):
        bad = self.good[:16] + b"\x00\x00\x00\x00" + self.good[20:]
        self.check_bad_header(bad, "Index offset 16824 but index length 0")

    def test_bad_index_eof(self):
        # Semi-random edit to the index offset value,
        bad = self.good[:13] + b"\x01" + self.good[14:]
        self.check_bad_header(
            bad,
            "Gap of 65536 bytes after final record end 16824, "
            "before 82360 where index starts?",
        )

    def test_bad_header_length(self):
        # 31 + 400 flow characters + 4 key bytes is 435, padded to 440
        self.assertEqual(self.good[24:26], b"\x01\xb8")
        for header_length in (432, 441, 448):
            bad = self.good[:24] + header_length.to_bytes(2, "big") + self.good[26:]
            self.check_bad_header(
                bad,
                f"SFF header length {header_length} does not match 31 + 400 flow "
                "characters + 4 key bytes padded to a multiple of 8",
            )

    def test_truncated_header(self):
        with self.assertRaises(ValueError) as cm:
            _sff_find_roche_index(BytesIO(self.good[:100]))
        self.assertEqual(
            str(cm.exception),
            "Premature end of file? Expected an SFF header of 440 bytes, "
            "but got only 100",
        )

    def test_truncated_header_parse(self):
        # 31 fixed bytes, 400 flow characters, 4 key bytes, 5 padding bytes;
        # cut inside the flow characters, the key, and the padding
        for size in (31, 100, 433, 437):
            with self.subTest(size=size):
                with warnings.catch_warnings():
                    # Before the fix, a misleading padding warning came first
                    warnings.simplefilter("error", BiopythonParserWarning)
                    self.check_bad_header(
                        self.good[:size],
                        "Premature end of file? Expected an SFF header of 440 "
                        f"bytes, but got only {size}",
                    )

    def test_negative_name_padding(self):
        # The first read's header says it is 16 bytes long, too short for
        # its 16 fixed bytes plus its 14 byte name
        self.assertEqual(self.good[440:444], b"\x00\x20\x00\x0e")
        bad = self.good[:440] + b"\x00\x10" + self.good[442:]
        msg = (
            "Malformed read header, says length is 16 but its fixed fields "
            "and 14 byte name need 30"
        )
        with self.assertRaises(ValueError) as cm:
            list(SeqIO.parse(BytesIO(bad), "sff"))
        self.assertEqual(str(cm.exception), msg)
        with self.assertRaises(ValueError) as cm:
            list(_sff_do_slow_index(BytesIO(bad)))
        self.assertEqual(str(cm.exception), msg)
        handle = BytesIO(bad)
        handle.seek(440)
        with self.assertRaises(ValueError) as cm:
            _sff_read_raw_record(handle, 400)
        self.assertEqual(str(cm.exception), msg)

    def test_write_records_twice(self):
        records = list(SeqIO.parse(BytesIO(self.good), "sff"))
        handle = BytesIO()
        writer = SffWriter(handle, index=False)
        writer.write_records(records[:2])
        written = handle.getvalue()
        with self.assertRaises(ValueError) as cm:
            writer.write_records(records[2:4])
        self.assertEqual(
            str(cm.exception),
            "SffWriter.write_records() can only be called once, "
            "as an SFF file has a single header",
        )
        self.assertEqual(handle.getvalue(), written)

    def test_slow_index_truncated_read_header(self):
        # File ends part way through the first read's name
        handle = BytesIO(self.good[: 440 + 16 + 5])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonParserWarning)
            with self.assertRaises(ValueError) as cm:
                list(_sff_do_slow_index(handle))
        self.assertEqual(
            str(cm.exception),
            "Premature end of SFF file, expected 14 byte read name for the "
            "read at offset 440 but got only 5",
        )

    def test_truncated_read(self):
        # The first read starts at 440: a 16 byte fixed header, a 14 byte
        # name, 2 bytes of padding, then a 1595 byte block of 400 flowgram
        # values (800 bytes) and 265 flow indexes, bases and qualities.
        # The last read starts at 15328 and its block is 1457 bytes.
        block = "flowgram and quality block"
        for size, offset, what, got in (
            (440, 440, "16 byte read header", 0),
            (445, 440, "16 byte read header", 5),
            (460, 440, "14 byte read name", 4),
            (480, 440, f"1595 byte {block}", 8),  # in the flowgram values
            (1900, 440, f"1595 byte {block}", 1428),  # in the qualities
            (16000, 15328, f"1457 byte {block}", 640),
        ):
            msg = (
                f"Premature end of SFF file, expected {what} for the read at "
                f"offset {offset} but got only {got}"
            )
            with self.subTest(size=size):
                with self.assertRaises(ValueError) as cm:
                    list(SeqIO.parse(BytesIO(self.good[:size]), "sff"))
                self.assertEqual(str(cm.exception), msg)
                handle = BytesIO(self.good[:size])
                handle.seek(offset)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", BiopythonParserWarning)
                    with self.assertRaises(ValueError) as cm:
                        _sff_read_raw_record(handle, 400)
                self.assertEqual(str(cm.exception), msg)

    def test_short_reads(self):
        # A raw or unbuffered handle may return fewer bytes than asked for
        # before the end of the file; that is not a truncated file.
        class ShortReads(BytesIO):
            def read(self, size=-1):
                if size is not None and size > 500:
                    size = 500
                return super().read(size)

        # The first read's 1595 byte flowgram and quality block is bigger
        # than any single read() on this handle returns.
        expected = next(SeqIO.parse(BytesIO(self.good), "sff"))
        record = next(SeqIO.parse(ShortReads(self.good), "sff"))
        self.assertEqual(record.id, expected.id)
        self.assertEqual(record.seq, expected.seq)
        self.assertEqual(record.letter_annotations, expected.letter_annotations)
        self.assertEqual(record.annotations, expected.annotations)
        raws = []
        for handle in (BytesIO(self.good), ShortReads(self.good)):
            handle.seek(440)
            raws.append(_sff_read_raw_record(handle, 400))
        self.assertEqual(raws[0], raws[1])

    def test_truncated_read_slow_index(self):
        block = "flowgram and quality block"
        for size, offset, what in (
            (440, 440, "16 byte read header"),
            (1900, 440, f"1595 byte {block}"),
            (16000, 15328, f"1457 byte {block}"),
        ):
            msg = (
                f"Premature end of SFF file, expected {what} for the read at "
                f"offset {offset}"
            )
            if size == 440:
                msg += " but got only 0"
            with self.subTest(size=size):
                with self.assertRaises(ValueError) as cm:
                    list(_sff_do_slow_index(BytesIO(self.good[:size])))
                self.assertEqual(str(cm.exception), msg)

    def test_write_len_mismatch(self):
        records = list(SeqIO.parse(BytesIO(self.good), "sff"))[:3]

        class Sized:
            # A sized iterable whose len() disagrees with what it yields
            def __len__(self):
                return 5

            def __iter__(self):
                return iter(records)

        with self.assertRaises(ValueError) as cm:
            SffWriter(BytesIO()).write_file(Sized())
        self.assertEqual(
            str(cm.exception),
            "Wrote 3 records, but len(records) was 5, so the read count "
            "in the SFF header is wrong",
        )

    def test_no_index(self):
        # Does a lot of work to create a no-index SFF file
        # (in the process checking this bit of SffWriter works)
        records = list(SeqIO.parse(BytesIO(self.good), "sff"))
        with BytesIO() as handle:
            writer = SffWriter(handle, index=False)
            count = writer.write_file(records)
            self.assertEqual(count, len(records))
            handle.seek(0)
            new = list(SeqIO.parse(handle, "sff"))
            self.assertEqual(len(records), len(new))
            for a, b in zip(records, new):
                self.assertEqual(a.id, b.id)
            handle.seek(0)
            with self.assertRaises(ValueError) as cm:
                values = _sff_find_roche_index(handle)
            err = str(cm.exception)
            self.assertEqual(err, "No index present in this SFF file")

    def test_unknown_index(self):
        # TODO - Add SFF file with no index,
        # self.assertEqual(str(err), "No index present in this SFF file")
        with open(
            support.DATA / "Roche" / "E3MFGYR02_alt_index_in_middle.sff", "rb"
        ) as handle:
            with self.assertRaises(ValueError) as cm:
                values = _sff_find_roche_index(handle)
        self.assertIn(
            str(cm.exception),
            (
                "Unknown magic number '.diy' in SFF index header:\n'.diy1.00'",
                "Unknown magic number b'.diy' in SFF index header:\nb'.diy1.00'",
            ),
        )

    def check_sff_read_roche_index(self, data, msg):
        handle = BytesIO(data)
        with self.assertRaises(ValueError) as cm:
            index = list(_sff_read_roche_index(handle))
        self.assertEqual(str(cm.exception), msg)

    def test_premature_end_of_index(self):
        self.check_sff_read_roche_index(self.good[:-50], "Premature end of file!")

    def test_index_name_no_null(self):
        self.assertEqual(self.good[17502:17503], b"\x00")
        self.check_sff_read_roche_index(
            self.good[:17502] + b"x" + self.good[17503:],
            "Expected a null terminator to the read name.",
        )

    def test_index_mft_version(self):
        self.assertEqual(self.good[16824:16832], b".mft1.00")
        self.check_sff_read_roche_index(
            self.good[:16828] + b"\x01\x02\x03\x04" + self.good[16832:],
            "Unsupported version in .mft index header, 1.2.3.4",
        )

    def test_index_mft_data_size(self):
        self.assertEqual(self.good[16824:16832], b".mft1.00")
        self.check_sff_read_roche_index(
            self.good[:16836] + b"\x00\x00\x00\x00" + self.good[16840:],
            "Problem understanding .mft index header, 764 != 8 + 8 + 548 + 0",
        )

    def test_index_lengths(self):
        # Reduce the number of reads from 10 to 9 so index loading fails...
        self.assertEqual(self.good[20:24], b"\x00\x00\x00\x0a")
        self.check_sff_read_roche_index(
            self.good[:20] + b"\x00\x00\x00\x09" + self.good[24:],
            "Problem with index length? 17568 vs 17588",
        )

    def test_no_manifest_xml(self):
        with open(support.DATA / "Roche" / "E3MFGYR02_no_manifest.sff", "rb") as handle:
            with self.assertRaises(ValueError) as cm:
                xml = ReadRocheXmlManifest(handle)
            self.assertEqual(str(cm.exception), "No XML manifest found")


class TestIndex(unittest.TestCase):
    """Test SFF index."""

    def test_manifest(self):
        filename = support.DATA / "Roche" / "E3MFGYR02_random_10_reads.sff"
        with open(filename, "rb") as handle:
            metadata = ReadRocheXmlManifest(handle)

    def test_both_ways(self):
        filename = support.DATA / "Roche" / "E3MFGYR02_random_10_reads.sff"
        with open(filename, "rb") as handle:
            index1 = sorted(_sff_read_roche_index(handle))
        with open(filename, "rb") as handle:
            index2 = sorted(_sff_do_slow_index(handle))
        self.assertEqual(index1, index2)
        with open(filename, "rb") as handle:
            self.assertEqual(len(index1), len(list(SffIterator(handle))))
        with open(filename, "rb") as handle:
            self.assertEqual(
                len(index1), len(list(SffIterator(BytesIO(handle.read()))))
            )


class TestAlternativeIndexes(unittest.TestCase):
    filename = support.DATA / "Roche" / "E3MFGYR02_random_10_reads.sff"
    with open(filename, "rb") as handle:
        sff = list(SffIterator(handle))

    def check_same(self, new_sff):
        self.assertEqual(len(self.sff), len(new_sff))
        for old, new in zip(self.sff, new_sff):
            self.assertEqual(old.id, new.id)
            self.assertEqual(old.seq, new.seq)

    def test_alt_index_at_end(self):
        with open(
            support.DATA / "Roche" / "E3MFGYR02_alt_index_at_end.sff", "rb"
        ) as handle:
            sff2 = list(SffIterator(handle))
        self.check_same(sff2)

    def test_alt_index_at_start(self):
        with open(
            support.DATA / "Roche" / "E3MFGYR02_alt_index_at_start.sff", "rb"
        ) as handle:
            sff2 = list(SffIterator(handle))
        self.check_same(sff2)

    def test_alt_index_in_middle(self):
        with open(
            support.DATA / "Roche" / "E3MFGYR02_alt_index_in_middle.sff", "rb"
        ) as handle:
            sff2 = list(SffIterator(handle))
        self.check_same(sff2)

    def test_index_at_start(self):
        with open(
            support.DATA / "Roche" / "E3MFGYR02_index_at_start.sff", "rb"
        ) as handle:
            sff2 = list(SffIterator(handle))
        self.check_same(sff2)

    def test_index_in_middle(self):
        with open(
            support.DATA / "Roche" / "E3MFGYR02_index_in_middle.sff", "rb"
        ) as handle:
            sff2 = list(SffIterator(handle))
        self.check_same(sff2)

    def test_trim(self):
        with open(self.filename, "rb") as handle:
            sff_trim = list(SffIterator(handle, trim=True))
        self.assertEqual(len(self.sff), len(sff_trim))
        for old, new in zip(self.sff, sff_trim):
            self.assertEqual(old.id, new.id)


class TestConcatenated(unittest.TestCase):
    """Test concatenated SFF files."""

    def test_parses_gzipped_stream(self):
        import gzip

        count = 0
        with gzip.open(
            support.DATA / "Roche" / "E3MFGYR02_random_10_reads.sff.gz", "rb"
        ) as fh:
            for record in SeqIO.parse(fh, "sff"):
                count += 1
        self.assertEqual(10, count)

    def test_parse1(self):
        count = 0
        caught = False
        try:
            for record in SeqIO.parse(
                support.DATA / "Roche" / "invalid_greek_E3MFGYR02.sff", "sff"
            ):
                count += 1
        except ValueError as err:
            self.assertIn(
                "Additional data at end of SFF file, perhaps "
                "multiple SFF files concatenated? "
                "See offset 65296",
                str(err),
                err,
            )
            caught = True
        self.assertTrue(caught, "Didn't spot concatenation")
        self.assertEqual(count, 24)

    def test_index1(self):
        with self.assertRaises(ValueError) as cm:
            d = SeqIO.index(
                support.DATA / "Roche" / "invalid_greek_E3MFGYR02.sff", "sff"
            )
        err = str(cm.exception)
        self.assertIn(
            "Additional data at end of SFF file, perhaps "
            "multiple SFF files concatenated? See offset 65296",
            err,
        )

    def test_parse2(self):
        count = 0
        caught = False
        try:
            for record in SeqIO.parse(
                support.DATA / "Roche" / "invalid_paired_E3MFGYR02.sff", "sff"
            ):
                count += 1
        except ValueError as err:
            self.assertIn(
                "Your SFF file is invalid, post index 5 byte "
                "null padding region ended '.sff' which could "
                "be the start of a concatenated SFF file? "
                "See offset 54371",
                str(err),
                err,
            )
            caught = True
        self.assertTrue(caught, "Didn't spot concatenation")
        self.assertEqual(count, 20)

    def test_index2(self):
        with self.assertRaises(ValueError) as cm:
            d = SeqIO.index(
                support.DATA / "Roche" / "invalid_paired_E3MFGYR02.sff", "sff"
            )
        self.assertIn(
            "Your SFF file is invalid, post index 5 byte "
            "null padding region ended '.sff' which could "
            "be the start of a concatenated SFF file? "
            "See offset 54371",
            str(cm.exception),
        )


class TestSelf(unittest.TestCase):
    """These tests were originally defined in SffIO.py as self-tests."""

    def test_read(self):
        filename = support.DATA / "Roche" / "E3MFGYR02_random_10_reads.sff"
        with open(filename, "rb") as handle:
            sff = list(SffIterator(handle))
        with open(filename, "rb") as handle:
            sff_trim = list(SffIterator(handle, trim=True))

        filename = support.DATA / "Roche" / "E3MFGYR02_random_10_reads_no_trim.fasta"
        fasta_no_trim = list(SeqIO.parse(filename, "fasta"))
        filename = support.DATA / "Roche" / "E3MFGYR02_random_10_reads_no_trim.qual"
        qual_no_trim = list(SeqIO.parse(filename, "qual"))

        filename = support.DATA / "Roche" / "E3MFGYR02_random_10_reads.fasta"
        fasta_trim = list(SeqIO.parse(filename, "fasta"))
        filename = support.DATA / "Roche" / "E3MFGYR02_random_10_reads.qual"
        qual_trim = list(SeqIO.parse(filename, "qual"))

        for s, sT, f, q, fT, qT in zip(
            sff, sff_trim, fasta_no_trim, qual_no_trim, fasta_trim, qual_trim
        ):
            self.assertEqual(len({s.id, f.id, q.id}), 1)  # All values are the same
            self.assertEqual(s.seq, f.seq)
            self.assertEqual(
                s.letter_annotations["phred_quality"],
                q.letter_annotations["phred_quality"],
            )
            self.assertEqual(
                len({s.id, sT.id, fT.id, qT.id}), 1
            )  # All values are the same
            self.assertEqual(sT.seq, fT.seq)
            self.assertEqual(
                sT.letter_annotations["phred_quality"],
                qT.letter_annotations["phred_quality"],
            )

    def test_write(self):
        filename = support.DATA / "Roche" / "E3MFGYR02_random_10_reads.sff"
        with open(filename, "rb") as handle:
            metadata = ReadRocheXmlManifest(handle)
        with open(filename, "rb") as handle:
            sff = list(SffIterator(handle))
        b_handle = BytesIO()
        w = SffWriter(b_handle, xml=metadata)
        w.write_file(sff)  # list
        data = b_handle.getvalue()
        # And again with an iterator...
        handle = BytesIO()
        w = SffWriter(handle, xml=metadata)
        w.write_file(iter(sff))
        self.assertEqual(data, handle.getvalue())
        # Check 100% identical to the original:
        with open(filename, "rb") as handle:
            original = handle.read()
        self.assertEqual(len(data), len(original))
        self.assertEqual(data, original)
        del data

    def test_index(self):
        filename = support.DATA / "Roche" / "greek.sff"
        with open(filename, "rb") as a_handle, open(filename, "rb") as b_handle:
            index1 = sorted(_sff_read_roche_index(a_handle))
            index2 = sorted(_sff_do_slow_index(b_handle))
            self.assertEqual(index1, index2)

    def test_read_wrong(self):
        filename = support.DATA / "Roche" / "greek.sff"
        with open(filename, "rb") as handle:
            self.assertRaises(ValueError, ReadRocheXmlManifest, handle)

        with open(filename, "rb") as handle:
            for record in SffIterator(handle):
                pass

            def fileiter(handle):
                for record in SffIterator(handle):
                    i = record.id

            self.assertRaises(ValueError, fileiter, handle)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)

if False:
    # Ugly code to make test files...
    index = ".diy1.00This is a fake index block (DIY = Do It Yourself), which is allowed under the SFF standard.\0"
    padding = len(index) % 8
    if padding:
        padding = 8 - padding
    index += chr(0) * padding
    assert len(index) % 8 == 0

    # Ugly bit of code to make a fake index at start
    index = ".diy1.00This is a fake index block (DIY = Do It Yourself), which is allowed under the SFF standard.\0"
    padding = len(index) % 8
    if padding:
        padding = 8 - padding
    index += chr(0) * padding
    with open(support.DATA / "Roche" / "E3MFGYR02_random_10_reads.sff", "rb") as handle:
        records = list(SffIterator(handle))
    with open(
        support.DATA / "Roche" / "E3MFGYR02_alt_index_at_start.sff", "w"
    ) as out_handle:
        w = SffWriter(out_handle, index=False, xml=None)
        # Fake the header...
        w._number_of_reads = len(records)
        w._index_start = 0
        w._index_length = 0
        w._key_sequence = records[0].annotations["flow_key"]
        w._flow_chars = records[0].annotations["flow_chars"]
        w._number_of_flows_per_read = len(w._flow_chars)
        w.write_header()
        w._index_start = out_handle.tell()
        w._index_length = len(index)
        out_handle.seek(0)
        w.write_header()  # this time with index info
        w.handle.write(index)
        for record in records:
            w.write_record(record)
    with open(
        support.DATA / "Roche" / "E3MFGYR02_alt_index_at_start.sff", "rb"
    ) as handle:
        records2 = list(SffIterator(handle))
    for old, new in zip(records, records2):
        assert old.seq == new.seq
    with open(
        support.DATA / "Roche" / "E3MFGYR02_alt_index_at_start.sff", "rb"
    ) as handle:
        i = list(_sff_do_slow_index(handle))

    # Ugly bit of code to make a fake index in middle
    index = ".diy1.00This is a fake index block (DIY = Do It Yourself), which is allowed under the SFF standard.\0"
    padding = len(index) % 8
    if padding:
        padding = 8 - padding
    index += chr(0) * padding
    with open(support.DATA / "Roche" / "E3MFGYR02_random_10_reads.sff", "rb") as handle:
        records = list(SffIterator(handle))
    with open(
        support.DATA / "Roche" / "E3MFGYR02_alt_index_in_middle.sff", "w"
    ) as out_handle:
        w = SffWriter(out_handle, index=False, xml=None)
        # Fake the header...
        w._number_of_reads = len(records)
        w._index_start = 0
        w._index_length = 0
        w._key_sequence = records[0].annotations["flow_key"]
        w._flow_chars = records[0].annotations["flow_chars"]
        w._number_of_flows_per_read = len(w._flow_chars)
        w.write_header()
        for record in records[:5]:
            w.write_record(record)
        w._index_start = out_handle.tell()
        w._index_length = len(index)
        w.handle.write(index)
        for record in records[5:]:
            w.write_record(record)
        out_handle.seek(0)
        w.write_header()  # this time with index info
    with open(
        support.DATA / "Roche" / "E3MFGYR02_alt_index_in_middle.sff", "rb"
    ) as handle:
        records2 = list(SffIterator(handle))
    for old, new in zip(records, records2):
        assert old.seq == new.seq
    with open(
        support.DATA / "Roche" / "E3MFGYR02_alt_index_in_middle.sff", "rb"
    ) as handle:
        j = list(_sff_do_slow_index(handle))

    # Ugly bit of code to make a fake index at end
    with open(support.DATA / "Roche" / "E3MFGYR02_random_10_reads.sff", "rb") as handle:
        records = list(SffIterator(handle))
    with open(
        support.DATA / "Roche" / "E3MFGYR02_alt_index_at_end.sff", "w"
    ) as out_handle:
        w = SffWriter(out_handle, index=False, xml=None)
        # Fake the header...
        w._number_of_reads = len(records)
        w._index_start = 0
        w._index_length = 0
        w._key_sequence = records[0].annotations["flow_key"]
        w._flow_chars = records[0].annotations["flow_chars"]
        w._number_of_flows_per_read = len(w._flow_chars)
        w.write_header()
        for record in records:
            w.write_record(record)
        w._index_start = out_handle.tell()
        w._index_length = len(index)
        out_handle.write(index)
        out_handle.seek(0)
        w.write_header()  # this time with index info
    with open(
        support.DATA / "Roche" / "E3MFGYR02_alt_index_at_end.sff", "rb"
    ) as handle:
        records2 = list(SffIterator(handle))
    for old, new in zip(records, records2):
        assert old.seq == new.seq
    with unittest.TestCase.assertRaises(None, ValueError):
        with open(
            support.DATA / "Roche" / "E3MFGYR02_alt_index_at_end.sff", "rb"
        ) as handle:
            print(ReadRocheXmlManifest(handle))
    with open(
        support.DATA / "Roche" / "E3MFGYR02_alt_index_at_end.sff", "rb"
    ) as handle:
        k = list(_sff_do_slow_index(handle))
    print("Done")
