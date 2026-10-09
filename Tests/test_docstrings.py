# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Run the docstring examples of every Bio and BioSQL module as doctests.

This file holds no tests itself. It anchors the doctest collector in
Tests/conftest.py, which imports each module by name and yields one test per
docstring, so ``pytest test_docstrings.py`` runs them all and
``pytest test_docstrings.py::Bio.Seq`` runs those of one module.
"""
