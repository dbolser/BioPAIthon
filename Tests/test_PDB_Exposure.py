# Copyright 2009-2011 by Eric Talevich.  All rights reserved.
# Revisions copyright 2009-2013 by Peter Cock.  All rights reserved.
# Revisions copyright 2013 Lenna X. Peterson. All rights reserved.
#
# Converted by Eric Talevich from an older unit test copyright 2002
# by Thomas Hamelryck.
#
# Merged related test files into one, by Joao Rodrigues (2020)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Unit tests for the Bio.PDB exposure classes."""

import unittest
import warnings

try:
    import numpy as np
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB."
    ) from None

import support

from Bio.PDB import CaPPBuilder
from Bio.PDB import ExposureCN
from Bio.PDB import HSExposureCA
from Bio.PDB import HSExposureCB
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning


class Exposure(unittest.TestCase):
    """Testing Bio.PDB.HSExposure."""

    def setUp(self):
        pdb_filename = support.DATA / "PDB" / "a_structure.pdb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            structure = PDBParser(PERMISSIVE=True).get_structure("X", pdb_filename)
        structure[1].detach_child("B")
        self.model = structure[1]
        # Look at first chain only
        a_residues = list(self.model["A"].child_list)
        self.assertEqual(86, len(a_residues))
        self.assertEqual(a_residues[0].get_resname(), "CYS")
        self.assertEqual(a_residues[1].get_resname(), "ARG")
        self.assertEqual(a_residues[2].get_resname(), "CYS")
        self.assertEqual(a_residues[3].get_resname(), "GLY")
        # ...
        self.assertEqual(a_residues[-3].get_resname(), "TYR")
        self.assertEqual(a_residues[-2].get_resname(), "ARG")
        self.assertEqual(a_residues[-1].get_resname(), "CYS")
        self.a_residues = a_residues
        self.radius = 13.0

    def test_HSExposureCA(self):
        """HSExposureCA."""
        _ = HSExposureCA(self.model, self.radius)
        residues = self.a_residues
        self.assertEqual(0, len(residues[0].xtra))
        self.assertEqual(0, len(residues[1].xtra))
        self.assertEqual(3, len(residues[2].xtra))
        self.assertAlmostEqual(
            0.81250973133184456, residues[2].xtra["EXP_CB_PCB_ANGLE"]
        )
        self.assertEqual(14, residues[2].xtra["EXP_HSE_A_D"])
        self.assertEqual(14, residues[2].xtra["EXP_HSE_A_U"])
        self.assertEqual(3, len(residues[3].xtra))
        self.assertAlmostEqual(1.3383737, residues[3].xtra["EXP_CB_PCB_ANGLE"])
        self.assertEqual(13, residues[3].xtra["EXP_HSE_A_D"])
        self.assertEqual(16, residues[3].xtra["EXP_HSE_A_U"])
        # ...
        self.assertEqual(3, len(residues[-2].xtra))
        self.assertAlmostEqual(
            0.77124014456278489, residues[-2].xtra["EXP_CB_PCB_ANGLE"]
        )
        self.assertEqual(24, residues[-2].xtra["EXP_HSE_A_D"])
        self.assertEqual(24, residues[-2].xtra["EXP_HSE_A_U"])
        self.assertEqual(0, len(residues[-1].xtra))

    def test_HSExposureCB(self):
        """HSExposureCB."""
        _ = HSExposureCB(self.model, self.radius)
        residues = self.a_residues
        self.assertEqual(0, len(residues[0].xtra))
        self.assertEqual(2, len(residues[1].xtra))
        self.assertEqual(20, residues[1].xtra["EXP_HSE_B_D"])
        self.assertEqual(5, residues[1].xtra["EXP_HSE_B_U"])
        self.assertEqual(2, len(residues[2].xtra))
        self.assertEqual(10, residues[2].xtra["EXP_HSE_B_D"])
        self.assertEqual(18, residues[2].xtra["EXP_HSE_B_U"])
        self.assertEqual(2, len(residues[3].xtra))
        self.assertEqual(7, residues[3].xtra["EXP_HSE_B_D"])
        self.assertEqual(22, residues[3].xtra["EXP_HSE_B_U"])
        # ...
        self.assertEqual(2, len(residues[-2].xtra))
        self.assertEqual(14, residues[-2].xtra["EXP_HSE_B_D"])
        self.assertEqual(34, residues[-2].xtra["EXP_HSE_B_U"])
        self.assertEqual(2, len(residues[-1].xtra))
        self.assertEqual(23, residues[-1].xtra["EXP_HSE_B_D"])
        self.assertEqual(15, residues[-1].xtra["EXP_HSE_B_U"])

    def test_HSExposureCB_glycine_missing_backbone_atom(self):
        """One unusable glycine must cost one residue, not the whole model."""
        # Record what an intact model gives, then start over.
        HSExposureCB(self.model, self.radius)
        expected = [dict(residue.xtra) for residue in self.a_residues]
        for residue in self.a_residues:
            residue.xtra.clear()

        # A glycine has no CB, so its pseudo-CB is built from N, C and CA.
        # Take the N away and it cannot be built at all.
        glycine = self.a_residues[3]
        self.assertEqual(glycine.get_resname(), "GLY")
        glycine.detach_child("N")

        HSExposureCB(self.model, self.radius)

        self.assertEqual(glycine.xtra, {})
        for i, residue in enumerate(self.a_residues):
            if residue is not glycine:
                self.assertEqual(residue.xtra, expected[i])

    def test_HSExposureCB_cb_on_ca(self):
        """A CB atom on its CA atom gives no half spheres for that residue."""
        HSExposureCB(self.model, self.radius)
        expected = [dict(residue.xtra) for residue in self.a_residues]
        for residue in self.a_residues:
            residue.xtra.clear()

        residue = self.a_residues[2]
        residue["CB"].coord = residue["CA"].coord.copy()
        HSExposureCB(self.model, self.radius)

        self.assertEqual(residue.xtra, {})
        for i, other in enumerate(self.a_residues):
            if other is not residue:
                self.assertEqual(other.xtra, expected[i])

    def test_HSExposureCA_cb_on_ca(self):
        """A CB atom on its CA atom gives no CB-pCB angle, but keeps the HSE."""
        HSExposureCA(self.model, self.radius)
        expected = [dict(residue.xtra) for residue in self.a_residues]
        for residue in self.a_residues:
            residue.xtra.clear()

        residue = self.a_residues[2]
        residue["CB"].coord = residue["CA"].coord.copy()
        HSExposureCA(self.model, self.radius)

        self.assertIsNone(residue.xtra["EXP_CB_PCB_ANGLE"])
        self.assertEqual(residue.xtra["EXP_HSE_A_U"], expected[2]["EXP_HSE_A_U"])
        self.assertEqual(residue.xtra["EXP_HSE_A_D"], expected[2]["EXP_HSE_A_D"])
        for i, other in enumerate(self.a_residues):
            if other is not residue:
                self.assertEqual(other.xtra, expected[i])

    def test_HSExposureCA_straight_cas(self):
        """Three CA atoms on a straight line give no pseudo-CB, so no HSE."""
        before, residue, after = self.a_residues[1:4]
        # Put the CA atoms on whole numbers of angstroms, so that the two
        # flanking CA atoms are exactly opposite and the pseudo-CB is zero.
        ca = residue["CA"].coord.round()
        residue["CA"].coord = ca
        before["CA"].coord = ca - [3, 0, 0]
        after["CA"].coord = ca + [3, 0, 0]
        peptides = CaPPBuilder().build_peptides(self.model)
        self.assertTrue(
            any(
                pp[i : i + 3] == [before, residue, after]
                for pp in peptides
                for i in range(len(pp))
            )
        )

        HSExposureCA(self.model, self.radius)

        self.assertEqual(residue.xtra, {})
        self.assertIn("EXP_HSE_A_U", self.a_residues[10].xtra)

    def test_HSExposureCB_coincident_ca(self):
        """A CA atom on top of this residue's CA is in neither half sphere."""
        residue = self.a_residues[2]
        # Move residues 60 and 61 together, so that they stay one peptide
        # for CaPPBuilder, until the CA of residue 60 is on that of residue 2.
        shift = residue["CA"].coord - self.a_residues[60]["CA"].coord
        for other in self.a_residues[60:62]:
            for atom in other:
                atom.coord = atom.coord + shift
        other = self.a_residues[60]
        other["CA"].coord = residue["CA"].coord.copy()  # exactly, not to roundoff
        self.assertEqual(other["CA"] - residue["CA"], 0)
        peptides = CaPPBuilder().build_peptides(self.model)
        self.assertTrue(any(other in pp for pp in peptides))

        HSExposureCB(self.model, self.radius)

        # Every other CA atom within the radius is counted once, up or down,
        # apart from the one on top of this residue's CA.
        neighbours = 0
        for pp in peptides:
            for ro in pp:
                if ro is not residue and 0 < ro["CA"] - residue["CA"] < self.radius:
                    neighbours += 1
        hse = residue.xtra["EXP_HSE_B_U"] + residue.xtra["EXP_HSE_B_D"]
        self.assertEqual(hse, neighbours)

    def test_ExposureCN(self):
        """HSExposureCN."""
        _ = ExposureCN(self.model, self.radius)
        residues = self.a_residues
        self.assertEqual(0, len(residues[0].xtra))
        self.assertEqual(1, len(residues[1].xtra))
        self.assertEqual(25, residues[1].xtra["EXP_CN"])
        self.assertEqual(1, len(residues[2].xtra))
        self.assertEqual(28, residues[2].xtra["EXP_CN"])
        self.assertEqual(1, len(residues[3].xtra))
        self.assertEqual(29, residues[3].xtra["EXP_CN"])
        # ...
        self.assertEqual(1, len(residues[-2].xtra))
        self.assertEqual(48, residues[-2].xtra["EXP_CN"])
        self.assertEqual(1, len(residues[-1].xtra))
        self.assertEqual(38, residues[-1].xtra["EXP_CN"])

    def test_negative_offset(self):
        """A negative or NaN offset raises ValueError, not assert."""
        for cls in (HSExposureCA, HSExposureCB, ExposureCN):
            for offset in (-1, float("nan")):
                with self.subTest(cls=cls.__name__, offset=offset):
                    with self.assertRaisesRegex(
                        ValueError, "offset must not be negative"
                    ):
                        cls(self.model, self.radius, offset=offset)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
