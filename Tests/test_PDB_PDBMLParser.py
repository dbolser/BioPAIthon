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

    def test_resolution(self):
        """The header gives the resolution recorded in the file."""
        pdbml_parser = PDBMLParser()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            for entry, resolution in [("1GBT", 2.0), ("3JQH", 2.201)]:
                with self.subTest(entry=entry):
                    structure = pdbml_parser.get_structure(
                        support.DATA / "PDB" / f"{entry}.xml"
                    )
                    self.assertEqual(structure.header["resolution"], resolution)

    def test_resolution_fallbacks(self):
        """Without refine's resolution, refine_hist's is read, then EM's."""
        with open(support.DATA / "PDB" / "3JQH.xml") as handle:
            data = handle.read()
        refine = "<PDBx:ls_d_res_high>2.201</PDBx:ls_d_res_high>"
        refine_hist = "<PDBx:d_res_high>2.201</PDBx:d_res_high>"
        category = "<PDBx:refine_histCategory>"
        for text in (refine, refine_hist, category):
            self.assertEqual(data.count(text), 1)
        nil = '<PDBx:ls_d_res_high xsi:nil="true" />'
        em = (
            "<PDBx:em_3d_reconstructionCategory>"
            '<PDBx:em_3d_reconstruction entry_id="3JQH" id="1">'
            "<PDBx:resolution>3.1</PDBx:resolution>"
            "</PDBx:em_3d_reconstruction>"
            "</PDBx:em_3d_reconstructionCategory>"
        )
        no_refine = data.replace(refine, "")
        neither = no_refine.replace(refine_hist, "")
        cases = [
            ("no refine", no_refine, 2.201),
            ("nil refine", data.replace(refine, nil), 2.201),
            ("em_3d_reconstruction", neither.replace(category, em + category), 3.1),
            ("none", neither, None),
        ]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            for case, text, resolution in cases:
                with self.subTest(case=case):
                    structure = PDBMLParser().get_structure(StringIO(text))
                    self.assertEqual(structure.header["resolution"], resolution)

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
