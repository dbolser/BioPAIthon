# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests for register_format in Bio.SeqIO and Bio.Align."""

import os
import subprocess
import sys
import tempfile
import types
import unittest
from io import StringIO

import support

from Bio import Align
from Bio import SeqIO
from Bio._io_registry import FormatRegistry
from Bio.Align import clustal
from Bio.SeqIO.FastaIO import FastaIterator
from Bio.SeqIO.FastaIO import FastaTwoLineIterator
from Bio.SeqIO.FastaIO import FastaTwoLineWriter
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.SeqIO.InsdcIO import GenBankIterator
from Bio.SeqIO.Interfaces import SequenceIterator
from Bio.SeqIO.SffIO import SffIterator
from Bio.SeqIO.SffIO import SffWriter
from Bio.SeqIO.UniprotIO import UniprotIterator

FASTA = support.DATA / "Fasta" / "f002"
FASTA_ONE = support.DATA / "Fasta" / "f001"
GENBANK = support.DATA / "GenBank" / "NC_005816.gb"
SFF = support.DATA / "Roche" / "E3MFGYR02_random_10_reads.sff"
UNIPROT = support.DATA / "SwissProt" / "multi_ex.xml"
CLUSTAL = support.DATA / "Clustalw" / "opuntia.aln"
BLAST_TAB = support.DATA / "Blast" / "tab_2226_tblastn_001.txt"


def restore(table, saved):
    """Put back a table's stored entries, as dict.copy saved them."""
    dict.clear(table)
    dict.update(table, saved)


class MarkingFastaIterator(FastaIterator):
    """FASTA parser which marks each record it makes."""

    def __next__(self):
        record = super().__next__()
        record.annotations["parsed by"] = "MarkingFastaIterator"
        return record


class MarkingGenBankIterator(GenBankIterator):
    """GenBank parser which marks each record it makes."""

    def __next__(self):
        record = super().__next__()
        record.annotations["parsed by"] = "MarkingGenBankIterator"
        return record


class NoMarkerIterator(SequenceIterator):
    """Text parser without record_start_marker."""

    modes = "t"

    def __next__(self):
        raise StopIteration


class BinaryIterator(SequenceIterator):
    """Binary parser with both indexing hooks."""

    modes = "b"
    record_start_marker = b">"

    @classmethod
    def parse_id_from_header(cls, line):
        return line[1:].split()[0].decode()

    def __next__(self):
        raise StopIteration


class NoHeaderRuleIterator(SequenceIterator):
    """Text parser with record_start_marker but the base parse_id_from_header."""

    modes = "t"
    record_start_marker = b">"

    def __next__(self):
        raise StopIteration


def parse_function(handle):
    """Parse FASTA, as a plain function."""
    return FastaIterator(handle)


class SeqIOTestCase(unittest.TestCase):
    """Restore Bio.SeqIO's tables after each test."""

    def setUp(self):
        for table in [
            SeqIO._FormatToIterator,
            SeqIO._FormatToWriter,
            SeqIO._converter,
        ]:
            self.addCleanup(restore, table, dict.copy(table))

    def write(self, records, fmt):
        handle = StringIO()
        SeqIO.write(records, handle, fmt)
        return handle.getvalue()

    def convert(self, filename, in_format, out_format):
        handle = StringIO()
        SeqIO.convert(filename, in_format, handle, out_format)
        return handle.getvalue()

    def index(self, filename, fmt):
        records = SeqIO.index(filename, fmt)
        self.addCleanup(records.close)
        return records

    def index_db(self, filename, fmt):
        records = SeqIO.index_db(":memory:", str(filename), fmt)
        self.addCleanup(records.close)
        return records


class SeqIORoundTrip(SeqIOTestCase):
    """A registered format works throughout Bio.SeqIO."""

    def check_round_trip(self, name):
        records = list(SeqIO.parse(FASTA, name))
        expected = list(SeqIO.parse(FASTA, "fasta"))
        self.assertEqual([r.id for r in records], [r.id for r in expected])
        self.assertEqual([r.seq for r in records], [r.seq for r in expected])
        record = SeqIO.read(FASTA_ONE, name)
        self.assertEqual(record.id, SeqIO.read(FASTA_ONE, "fasta").id)
        self.assertEqual(self.write(records, name), self.write(records, "fasta"))
        self.assertEqual(
            self.convert(GENBANK, "genbank", name),
            self.convert(GENBANK, "genbank", "fasta"),
        )
        self.assertEqual(
            self.convert(FASTA, name, "tab"), self.convert(FASTA, "fasta", "tab")
        )
        self.assertEqual(format(record, name), format(record, "fasta"))
        self.assertEqual(record.format(name), record.format("fasta"))

    def test_classes(self):
        SeqIO.register_format("test-fasta", FastaIterator, FastaWriter)
        self.check_round_trip("test-fasta")

    def test_spec_strings(self):
        SeqIO.register_format(
            "test-fasta",
            "Bio.SeqIO.FastaIO:FastaIterator",
            "Bio.SeqIO.FastaIO:FastaWriter",
        )
        self.check_round_trip("test-fasta")

    def test_roles_one_at_a_time(self):
        SeqIO.register_format("test-fasta", writer=FastaWriter)
        self.assertEqual(
            self.convert(FASTA, "fasta", "test-fasta"),
            self.convert(FASTA, "fasta", "fasta"),
        )
        with self.assertRaises(ValueError) as cm:
            SeqIO.parse(FASTA, "test-fasta")
        self.assertEqual(cm.exception.args, ("Unknown format 'test-fasta'",))
        SeqIO.register_format("test-fasta", iterator=FastaIterator)
        self.check_round_trip("test-fasta")


class SeqIOExistingNames(SeqIOTestCase):
    """A name's iterator or writer is replaced only with replace=True."""

    def test_builtin(self):
        for kwargs, role in [
            ({"iterator": MarkingFastaIterator}, "an iterator"),
            ({"writer": FastaTwoLineWriter}, "a writer"),
        ]:
            with self.subTest(role=role):
                with self.assertRaises(ValueError) as cm:
                    SeqIO.register_format("fasta", **kwargs)
                self.assertEqual(
                    str(cm.exception),
                    f"Format 'fasta' already has {role};"
                    " use replace=True to replace it",
                )
        self.assertIs(SeqIO._FormatToIterator["fasta"], FastaIterator)
        self.assertIs(SeqIO._FormatToWriter["fasta"], FastaWriter)

    def test_registered(self):
        SeqIO.register_format("test-fasta", FastaIterator)
        with self.assertRaises(ValueError):
            SeqIO.register_format("test-fasta", MarkingFastaIterator)
        self.assertIs(SeqIO._FormatToIterator["test-fasta"], FastaIterator)
        SeqIO.register_format("test-fasta", MarkingFastaIterator, replace=True)
        self.assertIs(SeqIO._FormatToIterator["test-fasta"], MarkingFastaIterator)

    def test_same_object_again(self):
        SeqIO.register_format("test-fasta", FastaIterator, FastaWriter)
        SeqIO.register_format("test-fasta", FastaIterator, FastaWriter)
        SeqIO.register_format("test-fasta", FastaIterator)
        self.assertIs(SeqIO._FormatToIterator["test-fasta"], FastaIterator)
        self.assertIs(SeqIO._FormatToWriter["test-fasta"], FastaWriter)

    def test_same_handler_as_spec_or_object(self):
        """A spec and the object it names count as the same, used or not."""
        spec = "Bio.SeqIO.FastaIO:FastaIterator"
        SeqIO.register_format("test-fasta", spec)
        SeqIO.register_format("test-fasta", FastaIterator)
        SeqIO.register_format("test-fasta", spec)
        self.assertEqual(dict.__getitem__(SeqIO._FormatToIterator, "test-fasta"), spec)
        self.assertIs(SeqIO._FormatToIterator["test-fasta"], FastaIterator)
        SeqIO.register_format("test-fasta", spec)
        SeqIO.register_format("test-fasta", FastaIterator)
        # The built-in table entries behave the same way:
        SeqIO.register_format("fasta", spec, "Bio.SeqIO.FastaIO:FastaWriter")
        SeqIO.register_format("fasta", FastaIterator, FastaWriter)
        with self.assertRaises(ValueError):
            SeqIO.register_format("test-fasta", "Bio.SeqIO.FastaIO:FastaWriter")

    def test_two_roles_all_or_nothing(self):
        SeqIO.register_format("test-fasta", FastaIterator)
        with self.assertRaises(ValueError):
            SeqIO.register_format("test-fasta", MarkingFastaIterator, FastaWriter)
        self.assertNotIn("test-fasta", SeqIO._FormatToWriter)
        SeqIO.register_format("test-other", writer=FastaWriter)
        with self.assertRaises(ValueError):
            SeqIO.register_format("test-other", FastaIterator, FastaTwoLineWriter)
        self.assertNotIn("test-other", SeqIO._FormatToIterator)
        self.assertIs(SeqIO._FormatToWriter["test-other"], FastaWriter)

    def test_bad_role_stores_nothing(self):
        for iterator, writer, error in [
            (FastaIterator, 42, TypeError),
            (42, FastaWriter, TypeError),
            (FastaIterator, "Bio.SeqIO.FastaIO", ValueError),
        ]:
            with self.subTest(iterator=iterator, writer=writer):
                with self.assertRaises(error):
                    SeqIO.register_format("test-fasta", iterator, writer)
                self.assertNotIn("test-fasta", SeqIO._FormatToIterator)
                self.assertNotIn("test-fasta", SeqIO._FormatToWriter)


class SeqIONames(SeqIOTestCase):
    """Names are checked exactly as SeqIO.parse checks them."""

    def test_same_errors_as_parse(self):
        for name in [None, 42, "", "FASTA", "Test-Fasta"]:
            with self.subTest(name=name):
                with self.assertRaises(Exception) as expected:
                    SeqIO.parse(StringIO(""), name)
                with self.assertRaises(Exception) as cm:
                    SeqIO.register_format(name, FastaIterator)
                self.assertIs(type(cm.exception), type(expected.exception))
                self.assertIn(type(cm.exception), (TypeError, ValueError))
                self.assertEqual(cm.exception.args, expected.exception.args)

    def test_no_roles(self):
        with self.assertRaises(TypeError):
            SeqIO.register_format("test-fasta")
        self.assertNotIn("test-fasta", SeqIO._FormatToIterator)
        self.assertNotIn("test-fasta", SeqIO._FormatToWriter)


class SeqIOReplace(SeqIOTestCase):
    """replace=True changes what parse, write, convert and index use."""

    def test_writer_used_by_write_convert_and_format(self):
        records = list(SeqIO.parse(GENBANK, "genbank"))
        unwrapped = self.write(records, "fasta-2line")
        self.assertNotEqual(self.write(records, "fasta"), unwrapped)
        SeqIO.register_format("fasta", writer=FastaTwoLineWriter, replace=True)
        self.assertEqual(self.write(records, "fasta"), unwrapped)
        self.assertEqual(self.convert(GENBANK, "genbank", "fasta"), unwrapped)
        self.assertEqual(self.convert(GENBANK, "gb", "fasta"), unwrapped)
        self.assertEqual(format(records[0], "fasta"), format(records[0], "fasta-2line"))
        # Only the shortcuts writing FASTA were dropped:
        self.assertNotIn(("genbank", "fasta"), SeqIO._converter)
        self.assertNotIn(("fastq", "fasta"), SeqIO._converter)
        self.assertIn(("fastq", "qual"), SeqIO._converter)

    def test_iterator_drops_shortcuts_reading_it(self):
        SeqIO.register_format("fastq", iterator=FastaIterator, replace=True)
        self.assertNotIn(("fastq", "fasta"), SeqIO._converter)
        self.assertNotIn(("fastq", "qual"), SeqIO._converter)
        self.assertIn(("fastq-sanger", "fasta"), SeqIO._converter)
        self.assertIn(("fastq-sanger", "fastq"), SeqIO._converter)

    def test_iterator_keeps_index_proxy(self):
        expected = [r.id for r in SeqIO.parse(FASTA, "fasta")]
        SeqIO.register_format("fasta", MarkingFastaIterator, replace=True)
        records = self.index(FASTA, "fasta")
        self.assertEqual(list(records), expected)
        for key in expected:
            self.assertEqual(
                records[key].annotations["parsed by"], "MarkingFastaIterator"
            )
        records = self.index_db(FASTA, "fasta")
        self.assertEqual(
            records[expected[0]].annotations["parsed by"], "MarkingFastaIterator"
        )

    def test_genbank_iterator_keeps_index_proxy(self):
        SeqIO.register_format("genbank", MarkingGenBankIterator, replace=True)
        records = self.index(GENBANK, "genbank")
        (key,) = records
        self.assertEqual(
            records[key].annotations["parsed by"], "MarkingGenBankIterator"
        )

    def test_iterator_without_record_start_marker(self):
        SeqIO.register_format("fasta", NoMarkerIterator, replace=True)
        with self.assertRaisesRegex(ValueError, "does not define record_start_marker"):
            SeqIO.index(FASTA, "fasta")

    def test_proxies_with_their_own_parser(self):
        """These index proxies refuse a replaced iterator, not a replaced writer."""
        for fmt, filename, builtin in [
            ("sff", SFF, SffIterator),
            ("sff-trim", SFF, SeqIO.SffIO._SffTrimIterator),
            ("uniprot-xml", UNIPROT, UniprotIterator),
        ]:
            with self.subTest(fmt=fmt):
                SeqIO.register_format(fmt, writer=FastaWriter, replace=True)
                self.assertTrue(len(self.index(filename, fmt)))
                replacement = type("Replacement", (builtin,), {})
                SeqIO.register_format(fmt, replacement, replace=True)
                with self.assertRaises(ValueError) as cm:
                    SeqIO.index(filename, fmt)
                self.assertEqual(str(cm.exception), f"Unsupported format {fmt!r}")
                with self.assertRaises(ValueError) as cm:
                    SeqIO.index_db(":memory:", str(filename), fmt)
                self.assertEqual(str(cm.exception), f"Unsupported format {fmt!r}")
                # Putting the built-in parser back makes it indexable again:
                SeqIO.register_format(fmt, builtin, replace=True)
                self.assertTrue(len(self.index(filename, fmt)))


class SeqIOIndexNewNames(SeqIOTestCase):
    """SeqIO.index and index_db scan a new format only if its parser can."""

    def test_qualifying_iterator(self):
        SeqIO.register_format("test-fasta", MarkingFastaIterator)
        expected = [r.id for r in SeqIO.parse(FASTA, "fasta")]
        records = self.index(FASTA, "test-fasta")
        self.assertEqual(list(records), expected)
        self.assertEqual(
            records[expected[1]].annotations["parsed by"], "MarkingFastaIterator"
        )
        records = self.index_db(FASTA, "test-fasta")
        self.assertEqual(sorted(records), sorted(expected))
        self.assertEqual(records[expected[2]].id, expected[2])

    def test_spec_string(self):
        SeqIO.register_format("test-fasta", "Bio.SeqIO.FastaIO:FastaIterator")
        self.assertEqual(len(self.index(FASTA, "test-fasta")), 3)

    def test_index_db_reload(self):
        """Reopening a database checks its format with the same rule."""
        SeqIO.register_format("test-fasta", FastaIterator)
        with tempfile.TemporaryDirectory() as directory:
            filename = os.path.join(directory, "test.idx")
            SeqIO.index_db(filename, str(FASTA), "test-fasta").close()
            records = SeqIO.index_db(filename)
            self.assertEqual(len(records), 3)
            records.close()
            SeqIO.register_format("test-fasta", parse_function, replace=True)
            with self.assertRaises(ValueError) as cm:
                SeqIO.index_db(filename)
            self.assertEqual(str(cm.exception), "Unsupported format 'test-fasta'")

    def test_unsupported(self):
        for iterator in [
            parse_function,
            BinaryIterator,
            NoHeaderRuleIterator,
            NoMarkerIterator,
        ]:
            with self.subTest(iterator=iterator):
                SeqIO.register_format("test-fasta", iterator, replace=True)
                with self.assertRaises(ValueError) as cm:
                    SeqIO.index(FASTA, "test-fasta")
                self.assertEqual(str(cm.exception), "Unsupported format 'test-fasta'")
                with self.assertRaises(ValueError) as cm:
                    SeqIO.index_db(":memory:", str(FASTA), "test-fasta")
                self.assertEqual(str(cm.exception), "Unsupported format 'test-fasta'")

    def test_unknown_name(self):
        with self.assertRaises(ValueError) as cm:
            SeqIO.index(FASTA, "test-absent")
        self.assertEqual(str(cm.exception), "Unsupported format 'test-absent'")

    def test_builtin_without_proxy_stays_unsupported(self):
        """A built-in name without an index proxy is never scanned for records."""
        hooks = {
            "record_start_marker": b">",
            "parse_id_from_header": vars(FastaIterator)["parse_id_from_header"],
        }
        for attribute, value in hooks.items():
            self.assertNotIn(attribute, vars(FastaTwoLineIterator))
            setattr(FastaTwoLineIterator, attribute, value)
            self.addCleanup(delattr, FastaTwoLineIterator, attribute)
        # The same class under a new name can be indexed:
        SeqIO.register_format("test-fasta", FastaTwoLineIterator)
        self.assertEqual(len(self.index(FASTA, "test-fasta")), 3)
        with self.assertRaises(ValueError) as cm:
            SeqIO.index(FASTA, "fasta-2line")
        self.assertEqual(str(cm.exception), "Unsupported format 'fasta-2line'")
        with self.assertRaises(ValueError) as cm:
            SeqIO.index_db(":memory:", str(FASTA), "fasta-2line")
        self.assertEqual(str(cm.exception), "Unsupported format 'fasta-2line'")


class AlignRegisterFormat(unittest.TestCase):
    """Bio.Align.register_format adds a case-insensitive format name."""

    def setUp(self):
        self.addCleanup(restore, Align._registry, dict.copy(Align._registry))
        self.formats = Align.formats
        self.alignment = Align.read(CLUSTAL, "clustal")
        self.expected = self.alignment.format("clustal")
        stream = StringIO()
        Align.write(self.alignment, stream, "clustal")
        self.expected_file = stream.getvalue()

    def tearDown(self):
        self.assertIs(Align.formats, self.formats)
        self.assertNotIn("star-aln", Align.formats)

    def check(self, *names):
        for name in names:
            with self.subTest(name=name):
                alignments = list(Align.parse(CLUSTAL, name))
                self.assertEqual(len(alignments), 1)
                self.assertEqual(alignments[0].format("clustal"), self.expected)
                alignment = Align.read(CLUSTAL, name)
                self.assertEqual(alignment.format("clustal"), self.expected)
                stream = StringIO()
                self.assertEqual(Align.write(alignment, stream, name), 1)
                self.assertEqual(stream.getvalue(), self.expected_file)
                self.assertEqual(alignment.format(name), self.expected)
                self.assertEqual(format(alignment, name), self.expected)

    def test_module(self):
        Align.register_format("STAR-ALN", clustal)
        self.assertIs(Align._registry["star-aln"], clustal)
        self.assertNotIn("STAR-ALN", Align._registry)
        self.check("star-aln", "Star-Aln", "STAR-ALN")

    def test_object_with_attributes(self):
        handler = types.SimpleNamespace(
            AlignmentIterator=clustal.AlignmentIterator,
            AlignmentWriter=clustal.AlignmentWriter,
        )
        Align.register_format("STAR-ALN", handler)
        self.check("star-aln", "Star-Aln")

    def test_spec_string(self):
        Align.register_format("STAR-ALN", "Bio.Align.clustal")
        self.check("star-aln", "Star-Aln")

    def test_reader_only(self):
        handler = types.SimpleNamespace(AlignmentIterator=clustal.AlignmentIterator)
        Align.register_format("star-aln", handler)
        self.assertEqual(len(list(Align.parse(CLUSTAL, "star-aln"))), 1)
        with self.assertRaises(ValueError):
            Align.write(self.alignment, StringIO(), "star-aln")

    def test_bad_arguments(self):
        for name, module, error in [
            (None, clustal, TypeError),
            (42, clustal, TypeError),
            ("", clustal, ValueError),
            ("star-aln", types.SimpleNamespace(), TypeError),
            ("star-aln", None, TypeError),
            ("FASTA", clustal, ValueError),
            ("fasta", clustal, ValueError),
        ]:
            with self.subTest(name=name, module=module):
                with self.assertRaises(error):
                    Align.register_format(name, module)
                self.assertNotIn("star-aln", Align._registry)
                self.assertEqual(list(Align._registry), list(Align.formats))

    def test_existing_names(self):
        Align.register_format("star-aln", clustal)
        Align.register_format("STAR-ALN", clustal)  # the same again: no-op
        handler = types.SimpleNamespace(AlignmentIterator=clustal.AlignmentIterator)
        with self.assertRaises(ValueError) as cm:
            Align.register_format("Star-Aln", handler)
        self.assertEqual(
            str(cm.exception),
            "Format 'star-aln' already exists; use replace=True to replace it",
        )
        Align.register_format("Star-Aln", handler, replace=True)
        self.assertIs(Align._registry["star-aln"], handler)
        Align.register_format("FASTA", clustal, replace=True)
        self.assertEqual(Align.read(CLUSTAL, "fasta").format("fasta"), self.expected)


class Registry(unittest.TestCase):
    """FormatRegistry.register and its built-in names."""

    def test_builtin(self):
        registry = FormatRegistry({"one": "a.b:c", "two": None})
        self.assertEqual(registry.builtin, frozenset(["one", "two"]))
        registry.register("three", "a.b:d")
        self.assertEqual(registry.builtin, frozenset(["one", "two"]))
        self.assertEqual(SeqIO._FormatToIterator.builtin, set(SeqIO._FormatToIterator))
        self.assertEqual(Align._registry.builtin, set(Align.formats))

    def test_register(self):
        marker = object()
        registry = FormatRegistry({"fmt": marker})
        registry.register("fmt", marker)
        registry.register("new", "a.b:c")
        registry.register("new", "a.b:c")
        for name, value in [("fmt", object()), ("new", "a.b:d"), ("fmt", None)]:
            with self.subTest(name=name, value=value):
                with self.assertRaises(ValueError) as cm:
                    registry.register(name, value)
                self.assertEqual(
                    str(cm.exception),
                    f"Format {name!r} already exists; use replace=True to replace it",
                )
        registry.register("fmt", None, replace=True)
        self.assertIsNone(dict.__getitem__(registry, "fmt"))
        self.assertEqual(dict.__getitem__(registry, "new"), "a.b:c")


class StaysWithoutImportlibMetadata(unittest.TestCase):
    """Using the built-in formats does not import importlib.metadata."""

    def test_builtin_formats(self):
        code = f"""\
import sys
from io import StringIO
from Bio import Align, AlignIO, SearchIO, SeqIO

records = list(SeqIO.parse({str(FASTA)!r}, "fasta"))
SeqIO.write(records, StringIO(), "fasta")
index = SeqIO.index({str(FASTA)!r}, "fasta")
index[records[0].id]
index.close()
db = SeqIO.index_db(":memory:", {str(FASTA)!r}, "fasta")
db[records[0].id]
db.close()
SeqIO.convert({str(FASTA)!r}, "fasta", StringIO(), "tab")
SeqIO.convert({str(GENBANK)!r}, "genbank", StringIO(), "fasta")
format(records[0], "fasta")
for fmt in ["clustal", "fasta"]:
    alignment = AlignIO.read({str(CLUSTAL)!r}, "clustal")
    AlignIO.write(alignment, StringIO(), fmt)
AlignIO.read(StringIO(format(alignment, "fasta")), "fasta")
Align.read({str(CLUSTAL)!r}, "clustal")
list(SearchIO.parse({str(BLAST_TAB)!r}, "blast-tab"))
assert "importlib.metadata" not in sys.modules, "importlib.metadata was imported"
"""
        result = subprocess.run(
            [sys.executable, "-W", "ignore", "-c", code],
            capture_output=True,
            text=True,
            cwd=support.DATA,
        )
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
