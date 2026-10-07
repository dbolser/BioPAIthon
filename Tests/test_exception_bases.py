# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests that exceptions for malformed input are ValueError subclasses.

Most of Biopython's parsers report bad input with ValueError or a subclass of
it (GenBank's ParserFailureError, SwissProt's SwissProtParserError and so on),
so ``except ValueError`` is how callers guard a parse.  The classes below used
to subclass Exception directly and slipped past that guard.
"""

import importlib
import unittest

# (module, class) for each exception that reports malformed input and was
# rebased from Exception onto ValueError.
INPUT_ERRORS = [
    ("Bio.CAPS", "AlignmentHasDifferentLengthsError"),
    ("Bio.Data.CodonTable", "TranslationError"),
    ("Bio.Nexus.Nexus", "NexusError"),
    ("Bio.Nexus.Trees", "TreeError"),
    ("Bio.PDB.PDBExceptions", "PDBConstructionException"),
    ("Bio.Phylo.NeXMLIO", "NeXMLError"),
    ("Bio.Phylo.NewickIO", "NewickError"),
    ("Bio.Phylo.PhyloXMLIO", "PhyloXMLError"),
]


class TestInputErrorBases(unittest.TestCase):
    """Each malformed-input exception can be caught as ValueError."""

    def test_input_errors_are_value_errors(self):
        for module_name, class_name in INPUT_ERRORS:
            with self.subTest(exception=f"{module_name}.{class_name}"):
                module = importlib.import_module(module_name)
                self.assertTrue(issubclass(getattr(module, class_name), ValueError))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
