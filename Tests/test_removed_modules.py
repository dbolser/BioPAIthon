# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests that removed modules explain themselves instead of vanishing.

``Bio.Application`` and the command line wrappers built on it were removed in
release 1.86.  Code that imports them is still common -- a GitHub code search
finds thousands of public files -- and without a stub each one gets a bare
``ModuleNotFoundError`` naming a module and nothing else.  ``Bio.Alphabet``,
removed in 1.78, has always been kept as a stub for this reason; these follow
it.

This is a deliberate difference from upstream Biopython, which ships no stub
for any of them.
"""

import unittest

# Each entry is a module and a phrase its error message must contain, so that
# a message cannot be reworded into uselessness without failing here.
REMOVED_MODULES = [
    ("Bio.Align.AlignInfo", "SummaryInfo"),
    ("Bio.Alphabet", "molecule_type"),
    ("Bio.Application", "subprocess"),
    ("Bio.Align.Applications", "MuscleCommandline"),
    ("Bio.Blast.Applications", "NcbiblastnCommandline"),
    ("Bio.Emboss.Applications", "WaterCommandline"),
    ("Bio.Phylo.Applications", "PhymlCommandline"),
    ("Bio.Sequencing.Applications", "Novoalign"),
    ("Bio.pairwise2", "PairwiseAligner"),
]


class TestRemovedModules(unittest.TestCase):
    """Importing a removed module raises ImportError, not ModuleNotFoundError."""

    def test_raises_import_error(self):
        for name, _ in REMOVED_MODULES:
            with self.subTest(module=name):
                with self.assertRaises(ImportError) as cm:
                    __import__(name)
                # A missing module raises ModuleNotFoundError, a subclass of
                # ImportError, so assertRaises alone would pass on no stub
                # at all.  The point of the stub is that it is *not* that.
                self.assertNotIsInstance(cm.exception, ModuleNotFoundError)

    def test_message_names_the_replacement(self):
        for name, phrase in REMOVED_MODULES:
            with self.subTest(module=name):
                with self.assertRaises(ImportError) as cm:
                    __import__(name)
                message = str(cm.exception)
                self.assertIn(name, message)
                self.assertIn(phrase, message)

    def test_still_importable_by_a_guarded_caller(self):
        """`except ImportError` around the old import must keep working."""
        for name, _ in REMOVED_MODULES:
            with self.subTest(module=name):
                try:
                    __import__(name)
                except ImportError:
                    continue
                self.fail(f"{name} imported without raising")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
