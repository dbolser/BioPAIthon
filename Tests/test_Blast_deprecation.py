# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for the import-time deprecation of Bio.Blast.NCBIWWW and NCBIXML.

Other test modules import these modules during collection, so an already
imported module would not warn again. Each test therefore drops the module
from sys.modules and imports it afresh, making the tests independent of
earlier imports in the same process.
"""

import importlib
import sys
import unittest
import warnings

from Bio import BiopythonDeprecationWarning


class DeprecationWarningTests(unittest.TestCase):
    """Each legacy module warns exactly once, naming its replacement."""

    def check_deprecation(self, module_name, replacement_fragments):
        """Freshly import module_name and return the warning message."""
        sys.modules.pop(module_name, None)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            importlib.import_module(module_name)
        deprecations = [
            item
            for item in caught
            if issubclass(item.category, BiopythonDeprecationWarning)
        ]
        self.assertEqual(
            len(deprecations),
            1,
            "importing %s should raise exactly one BiopythonDeprecationWarning,"
            " got %d" % (module_name, len(deprecations)),
        )
        message = str(deprecations[0].message)
        self.assertIn("%s has been deprecated" % module_name, message)
        for fragment in replacement_fragments:
            self.assertIn(fragment, message)

    def test_NCBIWWW(self):
        """Bio.Blast.NCBIWWW warns on import and points at Bio.Blast.qblast."""
        self.check_deprecation("Bio.Blast.NCBIWWW", ["Bio.Blast.qblast"])

    def test_NCBIXML(self):
        """Bio.Blast.NCBIXML warns on import and points at Bio.Blast read/parse."""
        self.check_deprecation(
            "Bio.Blast.NCBIXML", ["read and parse functions in Bio.Blast"]
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
