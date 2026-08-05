# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests that the re-exporting packages declare what they export.

These packages import names out of their submodules and offer them as their
own API -- ``from Bio.PDB import PDBParser`` is what the Tutorial teaches.
Without ``__all__`` that is an *implicit* re-export, which a downstream
``mypy --strict`` run rejects with "Module does not explicitly export
attribute".  Declaring ``__all__`` is what makes the promise explicit.
"""

import importlib
import inspect
import unittest

# Packages whose public API is re-exported from their submodules.  Note
# Bio.Restriction is deliberately absent: its ~1000 enzyme classes are
# synthesised at import and are invisible to a type checker either way, so
# it needs generated stubs rather than a list of names.
REEXPORTING_PACKAGES = [
    "Bio.Align",
    "Bio.Align.substitution_matrices",
    "Bio.AlignIO",
    "Bio.Blast",
    "Bio.CAPS",
    "Bio.Entrez",
    "Bio.ExPASy",
    "Bio.GenBank",
    "Bio.KEGG.Map",
    "Bio.PDB",
    "Bio.Pathway",
    "Bio.Phylo",
    "Bio.PopGen.GenePop",
    "Bio.SCOP",
    "Bio.SVDSuperimposer",
    "Bio.SearchIO.BlastIO",
    "Bio.SearchIO.ExonerateIO",
    "Bio.SearchIO.HHsuiteIO",
    "Bio.SearchIO.HmmerIO",
    "Bio.SearchIO.InfernalIO",
    "Bio.SearchIO.InterproscanIO",
    "Bio.SeqIO",
    "Bio.SeqUtils",
    "Bio.SwissProt",
    "Bio.TogoWS",
    "Bio.UniProt",
    "Bio.codonalign",
    "Bio.motifs",
    "Bio.motifs.jaspar",
    "Bio.phenotype",
]


class TestPublicExports(unittest.TestCase):
    """Each package declares __all__, and it matches what the package holds."""

    def test_all_is_declared(self):
        for name in REEXPORTING_PACKAGES:
            with self.subTest(package=name):
                module = importlib.import_module(name)
                self.assertIsInstance(getattr(module, "__all__", None), list)

    def test_every_exported_name_resolves(self):
        """A name in __all__ that is not there breaks `import *` at runtime."""
        for name in REEXPORTING_PACKAGES:
            module = importlib.import_module(name)
            for exported in module.__all__:
                with self.subTest(package=name, exported=exported):
                    self.assertTrue(hasattr(module, exported))

    def test_all_covers_the_public_namespace(self):
        """Nothing public is left undeclared.

        This is the guard that keeps the fix from decaying: adding a
        re-export without adding it to __all__ reintroduces exactly the
        implicit-re-export error this change removed.  Submodules are
        excluded because a type checker resolves those directly and does not
        require them to be declared.
        """
        for name in REEXPORTING_PACKAGES:
            module = importlib.import_module(name)
            public = {
                n
                for n in dir(module)
                if not n.startswith("_") and not inspect.ismodule(getattr(module, n))
            }
            with self.subTest(package=name):
                self.assertEqual(public, set(module.__all__))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
