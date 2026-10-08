# Copyright 2018 by Francesco Gastaldello. All rights reserved.
#
# Converted by Francesco Gastaldello from an older unit test copyright 2002
# by Thomas Hamelryck.
#
# Merged related test files into one, by Joao Rodrigues (2020)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Unit tests for the Bio.PDB.ResidueDepth module."""

import subprocess
import unittest
import warnings

import support

from Bio import BiopythonWarning
from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser
from Bio.PDB import ResidueDepth
from Bio.PDB.Atom import Atom
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.Residue import Residue
from Bio.PDB.ResidueDepth import _get_atom_radius


class MSMS_tests(unittest.TestCase):
    """Test calling MSMS via Bio.PDB.ResidueDepth."""

    @classmethod
    def setUpClass(cls):
        # Check if MSMS is installed
        try:
            v = subprocess.check_output(
                ["msms", "-h"], text=True, stderr=subprocess.STDOUT
            )
        except OSError:
            raise unittest.SkipTest(
                "Install MSMS if you want to use it from Biopython."
            )

        cls.pdbparser = PDBParser()
        cls.cifparser = MMCIFParser()

    def check_msms(self, prot_file, first_100_residues):
        """Wrap calls to MSMS and the respective tests."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            s = self.pdbparser.get_structure("X", prot_file)

        model = s[0]
        rd = ResidueDepth(model)
        residues = []
        for item in rd.property_list[:100]:
            residues.append(item[0].get_resname())
        self.assertEqual("".join(residues), first_100_residues)

    # def test_ResidueDepth_2XHE(self):
    #     self.check_msms(
    #         "PDB/2XHE.pdb",
    #         "HISMETSERLEULYSSERALAVALLYSTHRVALLEUTHRASNSERLEUARGSERVALALAASPGLYGLYASPTR"
    #         "PLYSVALLEUVALVALASPLYSPROALALEUARGMETILESERGLUCYSALAARGMETSERGLUILELEUASPL"
    #         "EUGLYVALTHRVALVALGLUASPVALSERLYSGLNARGLYSVALLEUPROGLNPHEHISGLYVALTYRPHEILE"
    #         "GLUPROTHRGLUGLUASNLEUASPTYRVALILEARGASPPHEALAASPARGTHRPROTHRTYRGLUALAALAHI"
    #         "SLEU",
    #     )

    def test_ResidueDepth_2BEG(self):
        self.check_msms(
            support.DATA / "PDB" / "2BEG.pdb",
            "LEUVALPHEPHEALAGLUASPVALGLYSERASNLYSGLYALAILEILEGLYLEUMETVALGLYGLYVALVALIL"
            "EALALEUVALPHEPHEALAGLUASPVALGLYSERASNLYSGLYALAILEILEGLYLEUMETVALGLYGLYVALV"
            "ALILEALALEUVALPHEPHEALAGLUASPVALGLYSERASNLYSGLYALAILEILEGLYLEUMETVALGLYGLY"
            "VALVALILEALALEUVALPHEPHEALAGLUASPVALGLYSERASNLYSGLYALAILEILEGLYLEUMETVALGL"
            "YGLY",
        )

    def test_ResidueDepth_1LCD(self):
        self.check_msms(
            support.DATA / "PDB" / "1LCD.pdb",
            "METLYSPROVALTHRLEUTYRASPVALALAGLUTYRALAGLYVALSERTYRGLNTHRVALSERARGVALVALAS"
            "NGLNALASERHISVALSERALALYSTHRARGGLULYSVALGLUALAALAMETALAGLULEUASNTYRILEPROA"
            "SNARG",
        )

    def test_ResidueDepth_1A8O(self):
        self.check_msms(
            support.DATA / "PDB" / "1A8O.pdb",
            "MSEASPILEARGGLNGLYPROLYSGLUPROPHEARGASPTYRVALASPARGPHETYRLYSTHRLEUARGALAGL"
            "UGLNALASERGLNGLUVALLYSASNTRPMSETHRGLUTHRLEULEUVALGLNASNALAASNPROASPCYSLYST"
            "HRILELEULYSALALEUGLYPROGLYALATHRLEUGLUGLUMSEMSETHRALACYSGLNGLY",
        )


class ResidueDepth_tests(unittest.TestCase):
    """Tests for Bio.PDB.ResidueDepth, except for running MSMS itself."""

    def test_pdb_to_xyzr(self):
        """Test generation of xyzr (atomic radii) file."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            p = PDBParser(PERMISSIVE=1)
            structure = p.get_structure("example", support.DATA / "PDB" / "1A8O.pdb")

        # Read radii produced with original shell script
        with open(support.DATA / "PDB" / "1A8O.xyzr") as handle:
            msms_radii = []
            for line in handle:
                fields = line.split()
                radius = float(fields[3])
                msms_radii.append(radius)

        model = structure[0]
        biopy_radii = []
        for atom in model.get_atoms():
            biopy_radii.append(_get_atom_radius(atom, rtype="united"))
        self.assertEqual(msms_radii, biopy_radii)

    @staticmethod
    def _make_atom(het, resname, name, element):
        residue = Residue((het, 1, " "), resname, "")
        atom = Atom(name, [0.0, 0.0, 0.0], 1.0, 1.0, " ", name, 1, element)
        residue.add(atom)
        return atom

    def test_atom_radius_rule_order(self):
        """Test radii where an earlier rule must win over a later one."""
        cases = [
            # het flag, residue, atom name, element, rtype, radius
            (" ", "ALA", "N", "N", "explicit", 1.54),
            (" ", "ALA", "N", "N", "united", 1.70),
            (" ", "ALA", "CA", "C", "united", 2.00),  # alpha carbon
            ("H_CA", "CA", "CA", "CA", "united", 1.97),  # calcium ion
            ("W", "HOH", "O", "O", "united", 1.60),  # water, not carbonyl O
            ("H_HEM", "HEM", "FE", "FE", "united", 1.30),  # heme iron
            ("H_FE2", "FE2", "FE", "FE", "united", 1.24),  # any other iron
        ]
        for het, resname, name, element, rtype, radius in cases:
            with self.subTest(resname=resname, name=name, rtype=rtype):
                atom = self._make_atom(het, resname, name, element)
                self.assertEqual(_get_atom_radius(atom, rtype=rtype), radius)

    def test_atom_radius_unknown(self):
        """Test the fallback for unknown atoms and the rtype check."""
        atom = self._make_atom("H_UNL", "UNL", "XX", "C")
        with self.assertWarns(BiopythonWarning):
            self.assertEqual(_get_atom_radius(atom), 0.01)
        with self.assertRaises(ValueError):
            _get_atom_radius(atom, rtype="bogus")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
