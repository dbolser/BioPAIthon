# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Offline tests for Bio.PDB.PDBList URL construction.

The download itself is faked, so these tests check which URL
``retrieve_pdb_file`` asks for, where it puts the result and what it
does when the download fails, without touching the network. The tests
that actually download from the PDB live in ``test_PDB_PDBList.py``.
"""

import contextlib
import gzip
import http.client
import importlib
import io
import os
import tempfile
import unittest
from unittest import mock

from Bio.PDB.PDBList import PDBList

pdb_list_module = importlib.import_module("Bio.PDB.PDBList")


class URLConstructionTests(unittest.TestCase):
    """Check the URLs retrieve_pdb_file builds and its failure path, mocked."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.pdblist = PDBList(
            pdb=self.tmp.name,
            obsolete_pdb=os.path.join(self.tmp.name, "obsolete"),
            verbose=False,
        )

    def retrieve(self, *args, **kwargs):
        """Call retrieve_pdb_file with urlretrieve faked; return (url, final path)."""
        requested = []

        def fake_urlretrieve(url, filename):
            requested.append(url)
            with gzip.open(filename, "wb") as handle:
                handle.write(b"payload")

        with mock.patch.object(pdb_list_module, "urlretrieve", fake_urlretrieve):
            with mock.patch.object(pdb_list_module, "urlcleanup"):
                final_file = self.pdblist.retrieve_pdb_file(*args, **kwargs)
        self.assertEqual(len(requested), 1)
        return requested[0], final_file

    def test_bcif_url(self):
        """BinaryCIF is fetched from the RCSB models server."""
        url, final_file = self.retrieve("127d", file_format="bcif")
        self.assertEqual(url, "https://models.rcsb.org/127d.bcif.gz")
        self.assertEqual(final_file, os.path.join(self.tmp.name, "27", "127d.bcif"))
        with open(final_file, "rb") as handle:
            self.assertEqual(handle.read(), b"payload")

    def test_mmcif_url(self):
        """mmCIF still comes from the configured wwPDB server."""
        url, final_file = self.retrieve("127d", file_format="mmCif")
        self.assertEqual(
            url,
            "https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/27/127d.cif.gz",
        )
        self.assertEqual(final_file, os.path.join(self.tmp.name, "27", "127d.cif"))

    def test_pdb_url(self):
        """PDB format still comes from the configured wwPDB server."""
        url, final_file = self.retrieve("127d", file_format="pdb")
        self.assertEqual(
            url,
            "https://files.wwpdb.org/pub/pdb/data/structures/divided/pdb/27/pdb127d.ent.gz",
        )
        self.assertEqual(final_file, os.path.join(self.tmp.name, "27", "pdb127d.ent"))

    def test_failed_download_returns_none(self):
        """A failed download returns None and prints why, leaving no file."""
        error = OSError("HTTP Error 404: Not Found")
        with mock.patch.object(pdb_list_module, "urlretrieve", side_effect=error):
            with mock.patch.object(pdb_list_module, "urlcleanup"):
                with contextlib.redirect_stdout(io.StringIO()) as output:
                    result = self.pdblist.retrieve_pdb_file("zzzz", file_format="pdb")
        self.assertIsNone(result)
        self.assertIn(
            "Desired structure not found or download failed. 'zzzz': "
            "HTTP Error 404: Not Found",
            output.getvalue(),
        )
        final_file = os.path.join(self.tmp.name, "zz", "pdbzzzz.ent")
        self.assertFalse(os.path.exists(final_file))

    def retrieve_failing(self, fake_urlretrieve, **kwargs):
        """Call retrieve_pdb_file for 127d in PDB format; return (result, output)."""
        with mock.patch.object(pdb_list_module, "urlretrieve", fake_urlretrieve):
            with mock.patch.object(pdb_list_module, "urlcleanup"):
                with contextlib.redirect_stdout(io.StringIO()) as output:
                    result = self.pdblist.retrieve_pdb_file(
                        "127d", file_format="pdb", **kwargs
                    )
        return result, output.getvalue()

    def test_non_gzip_response_leaves_no_file(self):
        """A non-gzip response returns None and leaves no file to reuse."""

        def fake_urlretrieve(url, filename):
            with open(filename, "wb") as handle:
                handle.write(b"<html>Please log in to the hotel wifi</html>")

        result, output = self.retrieve_failing(fake_urlretrieve)
        self.assertIsNone(result)
        self.assertIn(
            "Desired structure not found or download failed. '127d': ", output
        )
        self.assertEqual(os.listdir(os.path.join(self.tmp.name, "27")), [])
        # The next call must try again rather than report "Structure exists".
        url, final_file = self.retrieve("127d", file_format="pdb")
        with open(final_file, "rb") as handle:
            self.assertEqual(handle.read(), b"payload")

    def test_truncated_gzip_response_leaves_no_file(self):
        """A truncated gzip response returns None and leaves no file."""

        def fake_urlretrieve(url, filename):
            with open(filename, "wb") as handle:
                handle.write(gzip.compress(b"payload" * 100)[:-10])

        result, output = self.retrieve_failing(fake_urlretrieve)
        self.assertIsNone(result)
        self.assertIn("download failed. '127d'", output)
        self.assertEqual(os.listdir(os.path.join(self.tmp.name, "27")), [])

    def test_empty_response_leaves_no_file(self):
        """An empty response, or a gzip of nothing, returns None and no file."""
        for body in (b"", gzip.compress(b"")):
            with self.subTest(body=body):

                def fake_urlretrieve(url, filename, body=body):
                    with open(filename, "wb") as handle:
                        handle.write(body)

                result, output = self.retrieve_failing(fake_urlretrieve)
                self.assertIsNone(result)
                self.assertIn("download failed. '127d': empty download", output)
                self.assertEqual(os.listdir(os.path.join(self.tmp.name, "27")), [])

    def test_failed_overwrite_keeps_existing_file(self):
        """A failed download with overwrite=True leaves the old file intact."""
        url, final_file = self.retrieve("127d", file_format="pdb")

        def fake_urlretrieve(url, filename):
            with open(filename, "wb") as handle:
                handle.write(b"not gzip")

        result, output = self.retrieve_failing(fake_urlretrieve, overwrite=True)
        self.assertIsNone(result)
        self.assertEqual(os.listdir(os.path.join(self.tmp.name, "27")), ["pdb127d.ent"])
        with open(final_file, "rb") as handle:
            self.assertEqual(handle.read(), b"payload")

    def test_http_exception_returns_none(self):
        """A malformed server response (an HTTPException) returns None."""
        error = http.client.InvalidURL("nonnumeric port: 'x'")
        result, output = self.retrieve_failing(mock.Mock(side_effect=error))
        self.assertIsNone(result)
        self.assertIn(
            "Desired structure not found or download failed. '127d': "
            "nonnumeric port: 'x'",
            output,
        )
        self.assertFalse(
            os.path.exists(os.path.join(self.tmp.name, "27", "pdb127d.ent"))
        )

    def test_mmtf_raises(self):
        """The retired mmtf format is rejected before any download."""
        with self.assertRaises(ValueError) as context:
            self.pdblist.retrieve_pdb_file("127d", file_format="mmtf")
        self.assertIn("MMTF format is deprecated", str(context.exception))
        self.assertIn("bcif", str(context.exception))

    def test_obsolete_bcif_raises(self):
        """Obsolete structures cannot be requested in BinaryCIF format."""
        with self.assertRaises(ValueError) as context:
            self.pdblist.retrieve_pdb_file("347d", file_format="bcif", obsolete=True)
        self.assertEqual(
            "PDBList cannot retrieve obsolete structures in BinaryCIF format.",
            str(context.exception),
        )

    def test_invalid_format_raises(self):
        """An unknown format is rejected before any download."""
        with self.assertRaises(ValueError) as context:
            self.pdblist.retrieve_pdb_file("127d", file_format="mmtf2")
        self.assertIn("does not exist or is not supported", str(context.exception))

    def test_invalid_assembly_format_raises(self):
        """An unknown assembly format is rejected before any download."""
        with self.assertRaises(ValueError) as context:
            self.pdblist.retrieve_assembly_file("127d", 1, file_format="xml")
        self.assertIn("'xml' is not supported", str(context.exception))


class MalformedInputTests(unittest.TestCase):
    """Bad status files and directories raise ValueError, not assert."""

    def fake_urlopen(self, *lines):
        return mock.patch.object(
            pdb_list_module, "urlopen", return_value=io.BytesIO(b"".join(lines))
        )

    def test_status_list_bad_code(self):
        with self.fake_urlopen(b"1abc\n", b"12345\n"):
            with self.assertRaisesRegex(ValueError, "4 character PDB code"):
                PDBList.get_status_list("https://example.invalid/added.pdb")

    def test_obsolete_list_bad_code(self):
        with self.fake_urlopen(b"OBSLTE    14-JAN-04 1HKEX    1UUZ\n"):
            with self.assertRaisesRegex(ValueError, "4 character PDB code"):
                PDBList(verbose=False).get_all_obsolete()

    def test_update_pdb_without_local_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            missing = os.path.join(tmp, "missing")
            pdblist = PDBList(pdb=missing, verbose=False)
            with self.assertRaisesRegex(ValueError, "is not a directory"):
                pdblist.update_pdb()


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
