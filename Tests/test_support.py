# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for the shared test helpers in support.py."""

import sys
import unittest
import warnings

import support

import Bio.Blast
import Bio.Blast.NCBIXML
from Bio import BiopythonDeprecationWarning


class ImportFreshWithWarningsTests(unittest.TestCase):
    """import_fresh_with_warnings re-runs a module without replacing it."""

    def test_records_module_level_warnings(self):
        """The fresh import's warnings are returned, and it is a new object."""
        module, caught = support.import_fresh_with_warnings("Bio.Blast.NCBIXML")
        self.assertIsNot(module, Bio.Blast.NCBIXML)
        self.assertEqual(module.__name__, "Bio.Blast.NCBIXML")
        categories = [item.category for item in caught]
        self.assertIn(BiopythonDeprecationWarning, categories)

    def test_restores_module_and_parent_attribute(self):
        """sys.modules and the parent package hold the original object after."""
        old_module = sys.modules["Bio.Blast.NCBIXML"]
        self.assertIs(Bio.Blast.NCBIXML, old_module)
        support.import_fresh_with_warnings("Bio.Blast.NCBIXML")
        self.assertIs(sys.modules["Bio.Blast.NCBIXML"], old_module)
        self.assertIs(Bio.Blast.NCBIXML, old_module)

    def test_restores_on_import_failure(self):
        """State is restored even when the fresh import raises."""
        name = "Bio.Blast.NCBIXML"
        old_module = sys.modules[name]
        old_path = Bio.Blast.__path__
        Bio.Blast.__path__ = []
        try:
            with self.assertRaises(ModuleNotFoundError):
                support.import_fresh_with_warnings(name)
        finally:
            Bio.Blast.__path__ = old_path
        self.assertIs(sys.modules[name], old_module)
        self.assertIs(Bio.Blast.NCBIXML, old_module)

    def test_unimported_module_is_not_left_behind(self):
        """A module absent before the call is absent afterwards too."""
        name = "Bio.PopGen.GenePop.LargeFileParser"
        old_module = sys.modules.pop(name, None)
        import Bio.PopGen.GenePop

        if old_module is not None:
            delattr(Bio.PopGen.GenePop, "LargeFileParser")
        try:
            module, caught = support.import_fresh_with_warnings(name)
            self.assertEqual(module.__name__, name)
            self.assertNotIn(name, sys.modules)
            self.assertFalse(hasattr(Bio.PopGen.GenePop, "LargeFileParser"))
        finally:
            if old_module is not None:
                sys.modules[name] = old_module
                Bio.PopGen.GenePop.LargeFileParser = old_module

    def test_does_not_leak_warning_filters(self):
        """The "always" filter is confined to the call."""
        before = list(warnings.filters)
        support.import_fresh_with_warnings("Bio.Blast.NCBIXML")
        self.assertEqual(warnings.filters, before)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
