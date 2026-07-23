# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.File module."""

try:
    import sqlite3
except ImportError:
    # Run what tests we can in case sqlite3 was not installed
    sqlite3 = None

import os.path
import shutil
import struct
import tempfile
import unittest
from io import StringIO

from Bio import bgzf
from Bio import File
from Bio import MissingPythonDependencyError
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class RandomAccess(unittest.TestCase):
    """Random access tests."""

    def test_plain(self):
        """Test plain text file."""
        with File._open_for_random_access("Quality/example.fastq") as handle:
            self.assertIn("r", handle.mode)
            self.assertIn("b", handle.mode)

    def test_bgzf(self):
        """Test BGZF compressed file."""
        with File._open_for_random_access("Quality/example.fastq.bgz") as handle:
            self.assertIsInstance(handle, bgzf.BgzfReader)

    def test_gzip(self):
        """Test gzip compressed file."""
        self.assertRaises(
            ValueError, File._open_for_random_access, "Quality/example.fastq.gz"
        )


class AsHandleTestCase(unittest.TestCase):
    """Tests for as_handle function."""

    def setUp(self):
        """Initialise temporary directory."""
        # Create a directory to work in
        self.temp_dir = tempfile.mkdtemp(prefix="biopython-test")

    def tearDown(self):
        """Remove temporary directory."""
        shutil.rmtree(self.temp_dir)

    def _path(self, *args):
        return os.path.join(self.temp_dir, *args)

    def test_handle(self):
        """Test as_handle with a file-like object argument."""
        p = self._path("test_file.fasta")
        with open(p, "wb") as fp:
            with File.as_handle(fp) as handle:
                self.assertEqual(
                    fp,
                    handle,
                    "as_handle should return argument when given a file-like object",
                )
                self.assertFalse(handle.closed)

            self.assertFalse(
                handle.closed,
                "Exiting as_handle given a file-like object should not close the file",
            )

    def test_string_path(self):
        """Test as_handle with a string path argument."""
        p = self._path("test_file.fasta")
        mode = "wb"
        with File.as_handle(p, mode=mode) as handle:
            self.assertEqual(p, handle.name)
            self.assertEqual(mode, handle.mode)
            self.assertFalse(handle.closed)
        self.assertTrue(handle.closed)

    def test_path_object(self):
        """Test as_handle with a pathlib.Path object."""
        from pathlib import Path

        p = Path(self._path("test_file.fasta"))
        mode = "wb"
        with File.as_handle(p, mode=mode) as handle:
            self.assertEqual(str(p.absolute()), handle.name)
            self.assertEqual(mode, handle.mode)
            self.assertFalse(handle.closed)
        self.assertTrue(handle.closed)

    def test_custom_path_like_object(self):
        """Test as_handle with a custom path-like object."""

        class CustomPathLike:
            def __init__(self, path):
                self.path = path

            def __fspath__(self):
                return self.path

        p = CustomPathLike(self._path("test_file.fasta"))
        mode = "wb"
        with File.as_handle(p, mode=mode) as handle:
            self.assertEqual(p.path, handle.name)
            self.assertEqual(mode, handle.mode)
            self.assertFalse(handle.closed)
        self.assertTrue(handle.closed)

    def test_stringio(self):
        """Testing passing StringIO handles."""
        s = StringIO()
        with File.as_handle(s) as handle:
            self.assertIs(s, handle)


class BaseClassTests(unittest.TestCase):
    """Tests for _IndexedSeqFileProxy base class."""

    def test_instance_exception(self):
        self.assertRaises(TypeError, File._IndexedSeqFileProxy)

    def test_defaults_are_not_implemented(self):
        """The base class only promises to raise NotImplementedError."""

        class Concrete(File._IndexedSeqFileProxy):
            """Subclass which defers to the abstract base class throughout."""

            def __iter__(self):
                return super().__iter__()

            def get(self, offset):
                return super().get(offset)

        proxy = Concrete()
        self.assertRaises(NotImplementedError, iter, proxy)
        self.assertRaises(NotImplementedError, proxy.get, 0)
        # get_raw is documented as optional, and says so when unavailable
        with self.assertRaises(NotImplementedError) as cm:
            proxy.get_raw(0)
        self.assertIn("Not available for this file format", str(cm.exception))


class _StubHandle:
    """Stand-in for the binary handle a real proxy would hold open."""

    def __init__(self):
        self.closed = False

    def close(self):
        """Record that the handle was closed."""
        self.closed = True


class _StubProxy:
    """Minimal random access proxy driving _IndexedSeqFileDict directly.

    Real proxies re-read the sequence file, which makes it impossible to
    provoke a disagreement between the keys collected while indexing and the
    ids the parser later reports.  This stub separates the two, which is what
    happens in practice when an index outlives the file it describes.
    """

    def __init__(self, entries, records):
        self._handle = _StubHandle()
        self._entries = entries
        self._records = records

    def __iter__(self):
        """Return (key, offset, length) tuples as a real proxy would."""
        return iter(self._entries)

    def get(self, offset):
        """Return the SeqRecord registered for this offset."""
        return self._records[offset]


class IndexedSeqFileDictTests(unittest.TestCase):
    """Tests for the in-memory _IndexedSeqFileDict."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp(prefix="biopython-test")

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_str_of_empty_index(self):
        """An index over a file with no records prints as an empty dict."""
        path = os.path.join(self.temp_dir, "empty.fasta")
        with open(path, "w"):
            pass
        index = SeqIO.index(path, "fasta")
        try:
            self.assertEqual(len(index), 0)
            self.assertEqual(str(index), "{}")
        finally:
            index.close()

    def test_str_of_populated_index(self):
        """A populated index prints its first key and the object type."""
        index = SeqIO.index("Fasta/f002", "fasta")
        try:
            first_key = next(iter(index))
            self.assertEqual(str(index), "{%r : SeqRecord(...), ...}" % first_key)
        finally:
            index.close()

    def test_key_does_not_match_parsed_record(self):
        """Looking up a key the parsed record disagrees with is an error."""
        record = SeqRecord(Seq("ACGT"), id="Beta")
        proxy = _StubProxy([("Alpha", 0, 4)], {0: record})
        index = File._IndexedSeqFileDict(proxy, None, "repr", "SeqRecord")
        with self.assertRaises(ValueError) as cm:
            index["Alpha"]
        self.assertIn("Key did not match (Alpha vs Beta)", str(cm.exception))

    def test_close_closes_the_proxy_handle(self):
        """Closing the dictionary closes the handle the proxy holds."""
        record = SeqRecord(Seq("ACGT"), id="Alpha")
        proxy = _StubProxy([("Alpha", 0, 4)], {0: record})
        index = File._IndexedSeqFileDict(proxy, None, "repr", "SeqRecord")
        self.assertFalse(proxy._handle.closed)
        index.close()
        self.assertTrue(proxy._handle.closed)


def _proxy_factory(fmt, filename=None):
    """Return a proxy for the file, or whether the format is supported.

    This mirrors the closure Bio.SeqIO.index_db passes down to Bio.File.
    """
    from Bio.SeqIO._index import _FormatToRandomAccess

    if filename:
        return _FormatToRandomAccess[fmt](filename, fmt)
    return fmt in _FormatToRandomAccess


@unittest.skipIf(sqlite3 is None, "Requires sqlite3")
class SQLiteIndexTests(unittest.TestCase):
    """Tests for the SQLite backed _SQLiteManySeqFilesDict."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp(prefix="biopython-test")
        self.index_filename = os.path.join(self.temp_dir, "test.idx")

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def _tamper(self, key, value):
        """Overwrite one meta_data value in the index just built."""
        con = sqlite3.dbapi2.connect(self.index_filename)
        con.execute("UPDATE meta_data SET value = ? WHERE key = ?;", (value, key))
        con.commit()
        con.close()

    def test_missing_sqlite3_module(self):
        """Without the sqlite3 module index_db cannot work at all."""
        saved = File.sqlite3
        File.sqlite3 = None
        try:
            with self.assertRaises(MissingPythonDependencyError):
                SeqIO.index_db(self.index_filename, "Fasta/f002", "fasta")
        finally:
            File.sqlite3 = saved

    def test_unsupported_format_when_building(self):
        """Building an index for a format with no random access support fails."""
        with self.assertRaises(ValueError) as cm:
            SeqIO.index_db(self.index_filename, "Fasta/f002", "nonsense")
        self.assertIn("Unsupported format 'nonsense'", str(cm.exception))

    def test_repr(self):
        """The repr is the string handed in by Bio.SeqIO.index_db."""
        index = SeqIO.index_db(self.index_filename, "Fasta/f002", "fasta")
        try:
            self.assertEqual(
                repr(index),
                "SeqIO.index_db(%r, filenames=%r, format='fasta', key_function=None)"
                % (self.index_filename, ["Fasta/f002"]),
            )
        finally:
            index.close()

    def test_str_is_broken(self):
        """__str__ raises AttributeError; see the comment below.

        This is a bug, not intended behaviour.  _IndexedSeqFileDict.__str__
        formats using self._obj_repr, but _SQLiteManySeqFilesDict.__init__
        never sets that attribute, so printing an index_db dictionary blows
        up.  This test pins today's behaviour so that a fix has to update it
        deliberately; it is not an endorsement.
        """
        index = SeqIO.index_db(self.index_filename, "Fasta/f002", "fasta")
        try:
            with self.assertRaises(AttributeError) as cm:
                str(index)
            self.assertIn("_obj_repr", str(cm.exception))
        finally:
            index.close()

    def test_unfinished_database(self):
        """An index whose count is still -1 was never finished being built."""
        SeqIO.index_db(self.index_filename, "Fasta/f002", "fasta").close()
        self._tamper("count", -1)
        with self.assertRaises(ValueError) as cm:
            SeqIO.index_db(self.index_filename, "Fasta/f002", "fasta")
        self.assertIn("Unfinished/partial database", str(cm.exception))

    def test_corrupt_database_record_count(self):
        """The stored count must agree with the number of rows present."""
        index = SeqIO.index_db(self.index_filename, "Fasta/f002", "fasta")
        rows = len(index)
        index.close()
        self._tamper("count", rows + 7)
        with self.assertRaises(ValueError) as cm:
            SeqIO.index_db(self.index_filename, "Fasta/f002", "fasta")
        self.assertIn(
            "Corrupt database? %i entries not %i" % (rows, rows + 7), str(cm.exception)
        )

    def test_not_a_biopython_index(self):
        """An unrelated SQLite database is rejected, not misread."""
        con = sqlite3.dbapi2.connect(self.index_filename)
        con.execute("CREATE TABLE something_else (a TEXT);")
        con.commit()
        con.close()
        with self.assertRaises(ValueError) as cm:
            SeqIO.index_db(self.index_filename, "Fasta/f002", "fasta")
        self.assertIn("Not a Biopython index database?", str(cm.exception))

    def test_unsupported_format_when_reloading(self):
        """An index naming a format we cannot random access is rejected."""
        SeqIO.index_db(self.index_filename, "Fasta/f002", "fasta").close()
        self._tamper("format", "nonsense")
        with self.assertRaises(ValueError) as cm:
            # No format given, so it is taken from the index itself
            SeqIO.index_db(self.index_filename)
        self.assertIn("Unsupported format 'nonsense'", str(cm.exception))

    def test_reload_with_different_filenames(self):
        """Reloading an index with a different file list is rejected."""
        SeqIO.index_db(self.index_filename, "Fasta/f002", "fasta").close()
        with self.assertRaises(ValueError) as cm:
            SeqIO.index_db(self.index_filename, "Fasta/f001", "fasta")
        self.assertIn("Index file has different filenames", str(cm.exception))
        # New style indexes record paths relative to the index file, so the
        # message must not blame the old relative-to-$PWD behaviour.
        self.assertNotIn("original working directory", str(cm.exception))

    def test_get_raw_with_missing_key(self):
        """get_raw reports a missing key with KeyError, as documented."""
        index = SeqIO.index_db(self.index_filename, "Fasta/f002", "fasta")
        try:
            self.assertRaises(KeyError, index.get_raw, "no-such-key")
        finally:
            index.close()

    def test_stale_index_key_mismatch(self):
        """A key whose record has since changed identity is reported."""
        # Both versions of the file have identical record lengths, so the
        # offsets held in the index stay valid while the ids change.
        path = os.path.join(self.temp_dir, "seqs.fasta")
        with open(path, "w") as handle:
            handle.write(">Alpha\nACGT\n>Beta\nACGT\n")
        SeqIO.index_db(self.index_filename, path, "fasta").close()
        with open(path, "w") as handle:
            handle.write(">Gamma\nACGT\n>Delta\nACGT\n")
        index = SeqIO.index_db(self.index_filename, path, "fasta")
        try:
            with self.assertRaises(ValueError) as cm:
                index["Alpha"]
            self.assertIn("Key did not match (Alpha vs Gamma)", str(cm.exception))
        finally:
            index.close()


def open_descriptors(path):
    """Return how many of this process's descriptors point at the given file."""
    target = os.path.abspath(path)
    count = 0
    for name in os.listdir("/proc/self/fd"):
        try:
            if os.readlink(os.path.join("/proc/self/fd", name)) == target:
                count += 1
        except OSError:
            # Descriptor closed while we were looking at it
            pass
    return count


@unittest.skipIf(sqlite3 is None, "Requires sqlite3")
@unittest.skipUnless(os.path.isdir("/proc/self/fd"), "Requires /proc/self/fd")
class SQLiteIndexHandlePoolTests(unittest.TestCase):
    """Tests for the limit Bio.File puts on simultaneously open files.

    Bio.SeqIO.index_db does not expose max_open, so these build the
    dictionary directly with a small pool in order to reach the code which
    closes and reopens sequence files.
    """

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp(prefix="biopython-test")
        self.index_filename = os.path.join(self.temp_dir, "test.idx")

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def _index(self, filenames, fmt, max_open):
        return File._SQLiteManySeqFilesDict(
            self.index_filename,
            filenames,
            _proxy_factory,
            fmt,
            None,
            "repr",
            max_open=max_open,
        )

    def test_only_max_open_files_kept_open_while_building(self):
        """Building over more files than the pool allows closes the extras."""
        index = self._index(["Fasta/f001", "Fasta/f002"], "fasta", max_open=1)
        try:
            self.assertEqual(len(index._proxies), 1)
            self.assertEqual(len(index), 4)  # 1 record in f001, 3 in f002
            # The point of the pool is not to run out of file descriptors,
            # so the file dropped from it must really have been closed.
            self.assertEqual(open_descriptors("Fasta/f001"), 1)
            self.assertEqual(open_descriptors("Fasta/f002"), 0)
        finally:
            index.close()
        self.assertEqual(open_descriptors("Fasta/f001"), 0)

    def test_getitem_reopens_an_evicted_file(self):
        """Fetching from a closed file evicts a handle and reopens it."""
        index = self._index(["Fasta/f001", "Fasta/f002"], "fasta", max_open=1)
        try:
            keys = list(index)
            # Read from both files, forcing the pool to swap handles
            for key in keys:
                self.assertEqual(index[key].id, key)
            self.assertEqual(len(index._proxies), 1)
        finally:
            index.close()

    def test_get_raw_reopens_an_evicted_file(self):
        """get_raw can use the stored record length after reopening a file."""
        index = self._index(["Fasta/f001", "Fasta/f002"], "fasta", max_open=1)
        try:
            # Keys are ordered by file then offset, so the last one is the
            # last record of the second file: reading it needs both a fresh
            # handle and a seek to a non-zero offset.
            keys = list(index)
            raw = index.get_raw(keys[-1])
            self.assertTrue(raw.startswith(b">" + keys[-1].encode()))
            self.assertEqual(list(index._proxies), [1])
            self.assertEqual(open_descriptors("Fasta/f001"), 0)
            self.assertEqual(open_descriptors("Fasta/f002"), 1)
            # And back the other way
            raw = index.get_raw(keys[0])
            self.assertTrue(raw.startswith(b">" + keys[0].encode()))
            self.assertEqual(list(index._proxies), [0])
        finally:
            index.close()

    def test_get_raw_reopens_an_evicted_file_without_length(self):
        """SFF indexes hold no record length, so get_raw must ask the proxy."""
        # The Roche index in an SFF file gives offsets but no lengths, so the
        # length column is zero and Bio.File cannot take its usual shortcut.
        index = self._index(["Roche/greek.sff", "Roche/paired.sff"], "sff", max_open=1)
        try:
            # Keys come back ordered by file number, so the last one belongs
            # to the second file, which the pool will have closed.
            key = list(index)[-1]
            raw = index.get_raw(key)
            # An SFF read starts with a read header of two 16 bit lengths and
            # a 32 bit sequence length, then four 16 bit clip values, then the
            # read name.  The header length is padded to a multiple of eight.
            read_header_length, name_length, seq_len = struct.unpack(">2HI", raw[:8])
            self.assertEqual(read_header_length % 8, 0)
            self.assertEqual(raw[16 : 16 + name_length].decode(), key)
            self.assertEqual(len(raw) % 8, 0)
        finally:
            index.close()


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
