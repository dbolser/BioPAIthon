"""
Tests for the PDBML parser in the PDB package.

These tests rely on the principle that the structure returned by the PDBML parser should be the same as the structure
returned by the mmCIF parser for any PDB structure.
"""

import gzip
import unittest
import warnings
from io import StringIO

import support

from Bio.PDB import MMCIFParser
from Bio.PDB import PDBMLParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning


class TestPDBMLParser(unittest.TestCase):
    def test_get_structure(self):
        mmcif_parser = MMCIFParser()
        pdbml_parser = PDBMLParser()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            for entry in ["1GBT", "6WG6", "3JQH"]:
                mmcif_structure = mmcif_parser.get_structure(
                    entry, support.DATA / "PDB" / f"{entry}.cif"
                )
                if entry == "6WG6":
                    # This fixture is stored gzipped (39 MB uncompressed);
                    # the parser accepts a file handle, so read through gzip.
                    with gzip.open(support.DATA / "PDB" / f"{entry}.xml.gz") as source:
                        pdbml_structure = pdbml_parser.get_structure(source)
                else:
                    pdbml_structure = pdbml_parser.get_structure(
                        support.DATA / "PDB" / f"{entry}.xml"
                    )
                self.assertEqual(mmcif_structure, pdbml_structure)

    def test_get_structure_filehandle(self):
        mmcif_parser = MMCIFParser()
        pdbml_parser = PDBMLParser()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            for entry in ["1GBT"]:
                with (
                    open(support.DATA / "PDB" / f"{entry}.cif") as mmcif_file,
                    open(support.DATA / "PDB" / f"{entry}.xml") as pdbml_file,
                ):
                    mmcif_structure = mmcif_parser.get_structure(entry, mmcif_file)
                    pdbml_structure = pdbml_parser.get_structure(pdbml_file)
                self.assertEqual(mmcif_structure, pdbml_structure)

    def test_unexpected_atom_site_child(self):
        """A child of atom_siteCategory other than atom_site raises ValueError."""
        with open(support.DATA / "PDB" / "3JQH.xml") as handle:
            data = handle.read()
        old_start, old_end = '<PDBx:atom_site id="1">', "</PDBx:atom_site>"
        self.assertIn(old_start, data)
        data = data.replace(old_start, '<PDBx:other_site id="1">', 1)
        data = data.replace(old_end, "</PDBx:other_site>", 1)
        with self.assertRaisesRegex(
            ValueError,
            "Expected {http://pdbml.pdb.org/schema/pdbx-v50.xsd}atom_site in "
            "atom_siteCategory, found {http://pdbml.pdb.org/schema/pdbx-v50.xsd}"
            "other_site",
        ):
            PDBMLParser().get_structure(StringIO(data))


if __name__ == "__main__":
    unittest.main()
