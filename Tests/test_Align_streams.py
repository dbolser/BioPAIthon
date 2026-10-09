# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests for Bio.Align.parse reading streams, including ones that cannot seek.

A stream that cannot seek, such as a pipe or ``sys.stdin``, can be read once,
from start to end.  Streams that can seek are pinned too, as reading pipes
must not change how they behave.
"""

import contextlib
import io
import os
import tempfile
import threading
import unittest
from io import StringIO

import support

from Bio import Align

# One small fixture for each text format.  The big* formats are binary and
# need random access, so they are not included.
FIXTURES = {
    "a2m": support.DATA / "Clustalw" / "clustalw.a2m",
    "bed": support.DATA / "Blat" / "psl_34_001.bed",
    "chain": support.DATA / "Blat" / "psl_34_001.chain",
    "clustal": support.DATA / "Clustalw" / "clustalw.aln",
    "emboss": support.DATA / "Emboss" / "needle.txt",
    "exonerate": support.DATA / "Exonerate" / "exn_22_m_affine_local_cigar.exn",
    "fasta": support.DATA / "Clustalw" / "clustalw.fa",
    "hhr": support.DATA / "HHsuite" / "allx.hhr",
    "maf": support.DATA / "MAF" / "ucsc_test.maf",
    "mauve": support.DATA / "Mauve" / "combined.xmfa",
    "msf": support.DATA / "msf" / "W_prot.msf",
    "nexus": support.DATA / "Nexus" / "test_Nexus_input.nex",
    "phylip": support.DATA / "Phylip" / "horses.phy",
    "psl": support.DATA / "Blat" / "psl_34_001.psl",
    "sam": support.DATA / "Blat" / "psl_34_001.sam",
    "stockholm": support.DATA / "Stockholm" / "example.sth",
    "tabular": support.DATA / "Fasta" / "protein_m8CB.txt",
}

EMBOSS = FIXTURES["emboss"]  # five pairwise alignments


class ForwardOnlyHandle:
    """Mimic a network handle without seek and tell methods etc."""

    def __init__(self, handle):
        """Initialize the class."""
        self._handle = handle

    def __iter__(self):
        """Iterate."""
        return iter(self._handle)

    def __next__(self):
        """Get the next line."""
        return next(self._handle)

    def read(self, length=None):
        if length is None:
            return self._handle.read()
        else:
            return self._handle.read(length)

    def readline(self):
        return self._handle.readline()

    def close(self):
        return self._handle.close()


class WindowsStylePipe(ForwardOnlyHandle):
    """Mimic a pipe on Windows, which claims it can seek but cannot.

    On Windows, seekable() returns True for a pipe, and seek() silently does
    nothing (CPython issue gh-86768).  Wrapping a real pipe gives this class
    the pipe's file descriptor, as on Windows.
    """

    def seekable(self):
        return True

    def seek(self, offset, whence=0):
        return 0

    def tell(self):
        return 0

    def fileno(self):
        return self._handle.fileno()


@contextlib.contextmanager
def pipe_from(path):
    """Yield a text stream reading the contents of path through an OS pipe."""
    with open(path, "rb") as handle:
        data = handle.read()
    read_fd, write_fd = os.pipe()

    def feed():
        try:
            with open(write_fd, "wb") as stream:
                stream.write(data)
        except OSError:
            pass  # the reader closed its end before reading everything

    thread = threading.Thread(target=feed, daemon=True)
    thread.start()
    stream = open(read_fd)
    try:
        yield stream
    finally:
        stream.close()
        thread.join(timeout=10)


KINDS = ("forward-only", "pipe", "Windows-style pipe")


@contextlib.contextmanager
def nonseekable(kind, path):
    """Yield a text stream of the given kind, which cannot seek, reading path."""
    if kind == "forward-only":
        with open(path) as handle:
            yield ForwardOnlyHandle(handle)
    elif kind == "pipe":
        with pipe_from(path) as stream:
            yield stream
    elif kind == "Windows-style pipe":
        with pipe_from(path) as stream:
            yield WindowsStylePipe(stream)
    else:
        raise ValueError(kind)


def summarize(alignments):
    """Return the sequence ids and coordinates of each alignment.

    This iterates with ``for``, so given an alignments iterator it calls
    ``iter()``, which rewinds it.
    """
    return [
        (
            [sequence.id for sequence in alignment.sequences],
            alignment.coordinates.tolist(),
        )
        for alignment in alignments
    ]


class TestNonSeekableStreams(unittest.TestCase):
    """Parse every text format from streams that cannot seek."""

    def check(self, kind, use_list):
        for fmt, path in FIXTURES.items():
            with self.subTest(fmt=fmt):
                with Align.parse(path, fmt) as alignments:
                    expected = summarize(alignments)
                self.assertGreater(len(expected), 0)
                with nonseekable(kind, path) as stream:
                    alignments = Align.parse(stream, fmt)
                    if use_list:
                        alignments = list(alignments)
                    self.assertEqual(summarize(alignments), expected)

    def test_forward_only_for(self):
        self.check("forward-only", use_list=False)

    def test_forward_only_list(self):
        self.check("forward-only", use_list=True)

    def test_pipe_for(self):
        self.check("pipe", use_list=False)

    def test_pipe_list(self):
        self.check("pipe", use_list=True)

    def test_windows_style_pipe_for(self):
        self.check("Windows-style pipe", use_list=False)

    def test_windows_style_pipe_list(self):
        self.check("Windows-style pipe", use_list=True)


class TestNonSeekableBehaviour(unittest.TestCase):
    """Pin what a stream that cannot seek allows, and what it does not."""

    def test_empty(self):
        for kind in KINDS:
            for fmt in ("fasta", "a2m", "msf"):
                with self.subTest(kind=kind, fmt=fmt):
                    with nonseekable(kind, os.devnull) as stream:
                        alignments = Align.parse(stream, fmt)
                        with self.assertRaises(ValueError) as cm:
                            next(alignments)
                    self.assertEqual(str(cm.exception), "Empty file.")

    def test_len(self):
        for kind in KINDS:
            with self.subTest(kind=kind):
                with nonseekable(kind, EMBOSS) as stream:
                    alignments = Align.parse(stream, "emboss")
                    with self.assertRaises(TypeError):
                        len(alignments)
                    # len() read nothing, so all the alignments are still there,
                    self.assertEqual(len(list(alignments)), 5)
                    # and now that all of them have been read, len() knows:
                    self.assertEqual(len(alignments), 5)

    def test_second_pass(self):
        for kind in KINDS:
            with self.subTest(kind=kind):
                with nonseekable(kind, EMBOSS) as stream:
                    alignments = Align.parse(stream, "emboss")
                    self.assertEqual(len(list(alignments)), 5)
                    with self.assertRaises(io.UnsupportedOperation):
                        list(alignments)

    def test_for_after_next(self):
        for kind in KINDS:
            with self.subTest(kind=kind):
                with nonseekable(kind, EMBOSS) as stream:
                    alignments = Align.parse(stream, "emboss")
                    next(alignments)
                    with self.assertRaises(io.UnsupportedOperation):
                        for alignment in alignments:
                            pass

    def test_read(self):
        path = FIXTURES["phylip"]
        expected = summarize([Align.read(path, "phylip")])
        for kind in KINDS:
            with self.subTest(kind=kind):
                with nonseekable(kind, path) as stream:
                    alignment = Align.read(stream, "phylip")
                self.assertEqual(summarize([alignment]), expected)


class TestSeekableBehaviour(unittest.TestCase):
    """Pin how streams that can seek behave."""

    def test_empty(self):
        for fmt in ("fasta", "a2m", "msf"):
            with self.subTest(fmt=fmt):
                alignments = Align.parse(StringIO(""), fmt)
                with self.assertRaises(ValueError) as cm:
                    next(alignments)
                self.assertEqual(str(cm.exception), "Empty file.")

    def test_blank_lines(self):
        # A file of blank lines raises IndexError, not ValueError.  This is a
        # bug, pinned here only so that reading pipes does not change it.
        for fmt in ("fasta", "a2m", "msf"):
            with self.subTest(fmt=fmt):
                alignments = Align.parse(StringIO("\n\n"), fmt)
                with self.assertRaises(IndexError):
                    next(alignments)

    def test_at_end_of_file(self):
        fasta = ">a\nAC-GT\n>b\nACAGT\n"
        stream = StringIO(fasta)
        stream.seek(0, io.SEEK_END)
        alignments = Align.parse(stream, "fasta")
        with self.assertRaises(StopIteration):
            next(alignments)
        stream.seek(0, io.SEEK_END)
        alignments = Align.parse(stream, "fasta")
        self.assertEqual(
            summarize(alignments), [(["a", "b"], [[0, 2, 2, 4], [0, 2, 3, 5]])]
        )

    def test_for_after_next(self):
        with Align.parse(EMBOSS, "emboss") as alignments:
            first = summarize([next(alignments)])
            rest = summarize(alignments)  # a for loop, which starts again
        self.assertEqual(len(rest), 5)
        self.assertEqual(rest[:1], first)

    def test_second_pass(self):
        with open(EMBOSS) as handle:
            alignments = Align.parse(handle, "emboss")
            first = summarize(list(alignments))
            second = summarize(list(alignments))
            self.assertEqual(len(first), 5)
            self.assertEqual(second, first)
            self.assertEqual(len(alignments), 5)

    def test_fileno_raises_value_error(self):
        # A stream that can seek but has no working file descriptor.
        class Stream(StringIO):
            def fileno(self):
                raise ValueError("no file descriptor")

        with open(EMBOSS) as handle:
            stream = Stream(handle.read())
        alignments = Align.parse(stream, "emboss")
        self.assertEqual(len(alignments), 5)
        self.assertEqual(len(list(alignments)), 5)
        self.assertEqual(len(list(alignments)), 5)

    def test_spooled_temporary_file(self):
        # Its fileno() moves the data from memory to disk, and in text mode
        # loses it if a for loop over the stream has disabled tell().
        for fmt, path in FIXTURES.items():
            with self.subTest(fmt=fmt):
                with Align.parse(path, fmt) as alignments:
                    expected = summarize(alignments)
                with open(path) as handle:
                    data = handle.read()
                with tempfile.SpooledTemporaryFile(mode="w+") as stream:
                    stream.write(data)
                    stream.seek(0)
                    alignments = Align.parse(stream, fmt)
                    self.assertEqual(len(alignments), len(expected))
                    self.assertEqual(summarize(alignments), expected)
                    self.assertEqual(summarize(alignments), expected)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
