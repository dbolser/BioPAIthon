# Copyright 2017 by Bernhard Thiel.  All rights reserved.
# Revisions copyright 2024 James Krieger.

# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Unit tests for PARTS of the parse_pdb_header module of Bio.PDB."""

import unittest
from io import StringIO

try:
    import numpy as np  # noqa F401
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB."
    ) from None

import support

from Bio import BiopythonParserWarning
from Bio.PDB import PDBParser
from Bio.PDB.parse_pdb_header import _parse_remark_465
from Bio.PDB.parse_pdb_header import parse_pdb_header


class ParseReal(unittest.TestCase):
    """Testing with real PDB file(s)."""

    def test_1(self):
        """Parse the header of a known PDB file (1A8O)."""
        parser = PDBParser()
        struct = parser.get_structure("1A8O", support.DATA / "PDB" / "1A8O.pdb")
        self.assertAlmostEqual(struct.header["resolution"], 1.7)
        # Case-insensitive string comparisons
        known_strings = {
            "author": "T.R.Gamble,S.Yoo,F.F.Vajdos,U.K.Von Schwedler,D.K.Worthylake,H.Wang,J.P.Mccutcheon,W.I.Sundquist,C.P.Hill",
            "deposition_date": "1998-03-27",
            "head": "viral protein",
            "journal": "AUTH   T.R.GAMBLE,S.YOO,F.F.VAJDOS,U.K.VON SCHWEDLER,AUTH 2 D.K.WORTHYLAKE,H.WANG,J.P.MCCUTCHEON,W.I.SUNDQUIST,AUTH 3 C.P.HILLTITL   STRUCTURE OF THE CARBOXYL-TERMINAL DIMERIZATIONTITL 2 DOMAIN OF THE HIV-1 CAPSID PROTEIN.REF    SCIENCE                       V. 278   849 1997REFN                   ISSN 0036-8075PMID   9346481DOI    10.1126/SCIENCE.278.5339.849",
            "journal_reference": "t.r.gamble,s.yoo,f.f.vajdos,u.k.von schwedler, d.k.worthylake,h.wang,j.p.mccutcheon,w.i.sundquist, c.p.hill structure of the carboxyl-terminal dimerization domain of the hiv-1 capsid protein. science v. 278 849 1997 issn 0036-8075 9346481 10.1126/science.278.5339.849 ",
            "keywords": "capsid, core protein, hiv, c-terminal domain, viral protein",
            "name": "hiv capsid c-terminal domain",
            "release_date": "1998-10-14",
            "structure_method": "x-ray diffraction",
        }
        for key, expect in known_strings.items():
            self.assertEqual(struct.header[key].lower(), expect.lower())

    def test_2(self):
        """Parse the header of another PDB file (2BEG)."""
        parser = PDBParser()
        struct = parser.get_structure("2BEG", support.DATA / "PDB" / "2BEG.pdb")
        known_strings = {
            "author": "T.Luhrs,C.Ritter,M.Adrian,D.Riek-Loher,B.Bohrmann,H.Dobeli,D.Schubert,R.Riek",
            "deposition_date": "2005-10-24",
            "head": "protein fibril",
            "journal": "AUTH   T.LUHRS,C.RITTER,M.ADRIAN,D.RIEK-LOHER,B.BOHRMANN,AUTH 2 H.DOBELI,D.SCHUBERT,R.RIEKTITL   3D STRUCTURE OF ALZHEIMER'S AMYLOID-{BETA}(1-42)TITL 2 FIBRILS.REF    PROC.NATL.ACAD.SCI.USA        V. 102 17342 2005REFN                   ISSN 0027-8424PMID   16293696DOI    10.1073/PNAS.0506723102",
            "journal_reference": "t.luhrs,c.ritter,m.adrian,d.riek-loher,b.bohrmann, h.dobeli,d.schubert,r.riek 3d structure of alzheimer's amyloid-{beta}(1-42) fibrils. proc.natl.acad.sci.usa v. 102 17342 2005 issn 0027-8424 16293696 10.1073/pnas.0506723102 ",
            "keywords": "alzheimer's, fibril, protofilament, beta-sandwich, quenched hydrogen/deuterium exchange, pairwise mutagenesis, protein fibril",
            "name": "3d structure of alzheimer's abeta(1-42) fibrils",
            "release_date": "2005-11-22",
            "structure_method": "solution nmr",
        }
        for key, expect in known_strings.items():
            self.assertEqual(struct.header[key].lower(), expect.lower())

    def test_parse_pdb_with_remark_465(self):
        """Tests that parse_pdb_header now can identify some REMARK 465 entries."""
        header = parse_pdb_header(support.DATA / "PDB" / "2XHE.pdb")
        self.assertEqual(header["idcode"], "2XHE")
        self.assertTrue(header["has_missing_residues"])
        self.assertEqual(len(header["missing_residues"]), 142)
        self.assertIn(
            {
                "model": None,
                "res_name": "GLN",
                "chain": "B",
                "ssseq": 267,
                "insertion": None,
            },
            header["missing_residues"],
        )
        header = parse_pdb_header(support.DATA / "PDB" / "1A8O.pdb")
        self.assertFalse(header["has_missing_residues"])
        self.assertEqual(header["missing_residues"], [])

    def test_parse_remark_465(self):
        """A UNIT-test for the private function _parse_remark_465."""
        info = _parse_remark_465("GLU B   276")
        self.assertEqual(
            info,
            {
                "model": None,
                "res_name": "GLU",
                "chain": "B",
                "ssseq": 276,
                "insertion": None,
            },
        )

        info = _parse_remark_465("U R    -9")  # Based on PDB id: 7c7a
        self.assertEqual(
            info,
            {
                "model": None,
                "res_name": "U",
                "chain": "R",
                "ssseq": -9,
                "insertion": None,
            },
        )

        info = _parse_remark_465("2 GLU B   276B")
        self.assertEqual(
            info,
            {
                "model": 2,
                "res_name": "GLU",
                "chain": "B",
                "ssseq": 276,
                "insertion": "B",
            },
        )

        info = _parse_remark_465("A 2    11")
        self.assertEqual(
            info,
            {
                "model": None,
                "res_name": "A",
                "chain": "2",
                "ssseq": 11,
                "insertion": None,
            },
        )

        info = _parse_remark_465("1  DG B     9")
        self.assertEqual(
            info,
            {"model": 1, "res_name": "DG", "chain": "B", "ssseq": 9, "insertion": None},
        )

    def test_parse_header_line(self):
        """Unit test for parsing and converting fields in HEADER record."""
        header = parse_pdb_header(support.DATA / "PDB" / "header.pdb")
        self.assertEqual(header["head"], "structural genomics, unknown function")
        self.assertEqual(header["idcode"], "3EFG")
        self.assertEqual(header["deposition_date"], "2008-09-08")

    def test_parse_header_line_non_english_date(self):
        """Unit test for parsing and converting fields where date is not an english abbreviation."""
        with self.assertWarns(BiopythonParserWarning) as bad_date_warning:
            header = parse_pdb_header(
                support.DATA / "PDB" / "unrecognized_month_header.pdb", permissive=True
            )
        self.assertEqual(header["head"], "structural genomics, unknown function")
        self.assertEqual(header["idcode"], "3EFG")
        self.assertEqual(header["deposition_date"], "2008-00-08")
        self.assertEqual(
            str(bad_date_warning.warning),
            "Non-standard month in PDB header: Okt. Setting month to '00'.",
        )

    def test_parse_header_line_non_english_date_strict_mode(self):
        """Unit test for parsing and converting fields where date is not an english abbreviation."""
        with self.assertRaises(ValueError) as err:
            parse_pdb_header(
                support.DATA / "PDB" / "unrecognized_month_header.pdb", permissive=False
            )
        self.assertEqual(str(err.exception), "Non-standard month in PDB header: Okt.")

    def test_get_structure_non_english_date(self):
        """A non-english month in HEADER must not abort the whole structure parse."""
        with self.assertWarns(BiopythonParserWarning):
            structure = PDBParser().get_structure(
                "3EFG", support.DATA / "PDB" / "unrecognized_month_header.pdb"
            )
        self.assertEqual(structure.header["idcode"], "3EFG")
        self.assertEqual(structure.header["deposition_date"], "2008-00-08")
        self.assertEqual(len(list(structure.get_atoms())), 1)

    def test_get_structure_non_english_date_strict_mode(self):
        """With PERMISSIVE=False the error must name the offending month."""
        with self.assertRaises(ValueError) as err:
            PDBParser(PERMISSIVE=False).get_structure(
                "3EFG", support.DATA / "PDB" / "unrecognized_month_header.pdb"
            )
        self.assertEqual(str(err.exception), "Non-standard month in PDB header: Okt.")

    def test_parse_title_line(self):
        """Unit test for correct parsing of multiline title records."""
        header = parse_pdb_header(support.DATA / "PDB" / "1LCD.pdb")
        self.assertEqual(
            header["name"],
            "structure of the complex of lac repressor headpiece and an 11 "
            "base-pair half-operator determined by nuclear magnetic resonance "
            "spectroscopy and restrained molecular dynamics",
        )

    def test_parse_no_title(self):
        """Unit test for sensible result with no TITLE line."""
        header = parse_pdb_header(support.DATA / "PDB" / "occupancy.pdb")
        self.assertEqual(header["name"], "")

    def test_compound_chain_ids_keep_case(self):
        """COMPND CHAIN values keep their case, as "A" and "a" are different chains."""
        # The COMPND block of PDB entry 7DK4 (biopython/biopython#4176), with
        # a wrapped CHAIN list added to molecule 3, then two atoms.
        lines = [
            "HEADER    VIRAL PROTEIN                           23-NOV-20   7DK4              \n",
            "COMPND    MOL_ID: 1;                                                            \n",
            "COMPND   2 MOLECULE: THE HEAVY CHAIN OF 2H2 FAB;                                \n",
            "COMPND   3 CHAIN: a, c, e;                                                      \n",
            "COMPND   4 MOL_ID: 2;                                                           \n",
            "COMPND   5 MOLECULE: THE LIGHT CHAIN OF 2H2 FAB;                                \n",
            "COMPND   6 CHAIN: b, d, f;                                                      \n",
            "COMPND   7 MOL_ID: 3;                                                           \n",
            "COMPND   8 MOLECULE: SPIKE GLYCOPROTEIN;                                        \n",
            "COMPND   9 CHAIN: A, B, C,                                                      \n",
            "COMPND  10 D, g;                                                                \n",
            "COMPND  11 SYNONYM: S GLYCOPROTEIN,E2,PEPLOMER PROTEIN,SARS-COV-2 S PROTEIN;    \n",
            "COMPND  12 ENGINEERED: YES                                                      \n",
            "ATOM      1  N   GLN A  14     215.703 239.017 162.981  1.00 87.37           N  \n",
            "ATOM      2  N   GLN a   1     190.432 205.226 241.016  1.00 99.94           N  \n",
            "END                                                                             \n",
        ]
        header = parse_pdb_header(StringIO("".join(lines)))
        compound = header["compound"]
        # MOL_ID dispatch and the lowercase key names are unchanged
        self.assertEqual(list(compound), ["1", "2", "3"])
        self.assertEqual(
            list(compound["3"]),
            ["misc", "molecule", "chain", "synonym", "engineered"],
        )
        # The chain IDs are exactly as written in the file
        self.assertEqual(compound["1"]["chain"], "a, c, e")
        self.assertEqual(compound["2"]["chain"], "b, d, f")
        self.assertEqual(
            [c.strip() for c in compound["3"]["chain"].split(",")],
            ["A", "B", "C", "D", "g"],
        )
        # Free text is still lowercased, as before
        self.assertEqual(compound["3"]["molecule"], "spike glycoprotein")
        self.assertEqual(compound["3"]["engineered"], "yes")
        # So each chain of the structure finds its own molecule
        structure = PDBParser().get_structure("7DK4", StringIO("".join(lines)))
        chain_to_mol = {
            chain_id.strip(): mol_id
            for mol_id, info in structure.header["compound"].items()
            for chain_id in info.get("chain", "").split(",")
            if chain_id.strip()
        }
        self.assertEqual(
            {chain.id: chain_to_mol[chain.id] for chain in structure[0]},
            {"A": "3", "a": "1"},
        )

    def test_parse_pdb_with_remark_99(self):
        """Tests that parse_pdb_header can identify REMARK 99 ASTRAL entries."""
        header = parse_pdb_header(support.DATA / "PDB" / "d256ba_.ent")
        self.assertIn("astral", header)
        self.assertEqual(header["astral"]["SCOP-sccs"], "a.24.3.1")
        self.assertEqual(header["astral"]["Source-PDB"], "256b")
        self.assertEqual(header["astral"]["Region"], "a:")
        self.assertEqual(header["astral"]["ASTRAL-SPACI"], "0.72")

    def test_parse_pdb_with_remark_350_biomoltrans(self):
        """Tests that parse_pdb_header now can identify some REMARK 350 entries."""
        header = parse_pdb_header(support.DATA / "PDB" / "2XHE.pdb")
        self.assertEqual(
            header["biomoltrans"],
            {
                "1": [
                    ["A", "B"],
                    "  1.000000  0.000000  0.000000        0.00000            \n",
                    "  0.000000  1.000000  0.000000        0.00000            \n",
                    "  0.000000  0.000000  1.000000        0.00000            \n",
                ]
            },
        )

    def test_parse_pdb_without_remark_350_biomoltrans(self):
        header = parse_pdb_header(support.DATA / "PDB" / "1LCD.pdb")
        self.assertEqual(header["biomoltrans"], {})


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
