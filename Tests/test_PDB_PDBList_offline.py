# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Offline tests for Bio.PDB.PDBList URL construction.

The download itself is faked, so these tests check which URL
``retrieve_pdb_file`` asks for and where it puts the result, without
touching the network. The tests that actually download from the PDB
live in ``test_PDB_PDBList.py``.
"""

import gzip
import os
import tempfile
import unittest
from unittest import mock

from Bio.PDB.PDBList import PDBList


class URLConstructionTests(unittest.TestCase):
    """Check the URLs retrieve_pdb_file builds, with the download mocked."""

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

        with mock.patch("Bio.PDB.PDBList.urlretrieve", fake_urlretrieve):
            with mock.patch("Bio.PDB.PDBList.urlcleanup"):
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


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
