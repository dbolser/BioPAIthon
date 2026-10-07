# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO fasta-m10 indexing."""

import unittest

import support
from search_tests_common import CheckIndex


class FastaM10IndexCases(CheckIndex):
    fmt = "fasta-m10"

    def test_output_002(self):
        """Test fasta-m10 indexing, fasta34, multiple queries."""
        filename = support.DATA / "Fasta" / "output002.m10"
        self.check_index(filename, self.fmt)

    def test_output_001(self):
        """Test fasta-m10 indexing, fasta35, multiple queries."""
        filename = support.DATA / "Fasta" / "output001.m10"
        self.check_index(filename, self.fmt)

    def test_output_005(self):
        """Test fasta-m10 indexing, ssearch35, multiple queries."""
        filename = support.DATA / "Fasta" / "output005.m10"
        self.check_index(filename, self.fmt)

    def test_output_008(self):
        """Test fasta-m10 indexing, tfastx36, multiple queries."""
        filename = support.DATA / "Fasta" / "output008.m10"
        self.check_index(filename, self.fmt)

    def test_output_009(self):
        """Test fasta-m10 indexing, fasta36, multiple queries."""
        filename = support.DATA / "Fasta" / "output009.m10"
        self.check_index(filename, self.fmt)

    def test_output_010(self):
        """Test fasta-m10 indexing, fasta36, single query, no hits."""
        filename = support.DATA / "Fasta" / "output010.m10"
        self.check_index(filename, self.fmt)

    def test_output_011(self):
        """Test fasta-m10 indexing, fasta36, single query, hits with single hsp."""
        filename = support.DATA / "Fasta" / "output011.m10"
        self.check_index(filename, self.fmt)

    def test_output_012(self):
        """Test fasta-m10 indexing, fasta36, single query with multiple hsps."""
        filename = support.DATA / "Fasta" / "output012.m10"
        self.check_index(filename, self.fmt)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
