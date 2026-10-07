# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Offline tests for Bio.UniProt search results.

The REST responses are faked, so these tests run without the network. The
tests that query UniProt itself live in ``test_UniProt.py``.
"""

import io
import json
import unittest
from unittest import mock

from Bio import UniProt


class FakeResponse(io.BytesIO):
    """A UniProt search response holding one batch and no next link."""

    def __init__(self, results, total):
        """Hold the JSON for results, and a header giving the total."""
        super().__init__(json.dumps({"results": results}).encode())
        self.headers = {"x-total-results": str(total)}


def fake_search(results, total):
    """Return search results built from one faked response."""
    response = FakeResponse(results, total)
    with mock.patch.object(UniProt, "urlopen", return_value=response):
        return UniProt._UniProtSearchResults("https://rest.uniprot.org/fake")


class SearchResultsTests(unittest.TestCase):
    def test_fewer_results_than_reported(self):
        """Running out of batches early raises ValueError, not assert."""
        results = fake_search([{"id": 1}, {"id": 2}], total=3)
        iterator = iter(results)
        self.assertEqual(next(iterator), {"id": 1})
        self.assertEqual(next(iterator), {"id": 2})
        with self.assertRaisesRegex(ValueError, "reported 3 results, but gave no"):
            next(iterator)

    def test_slices(self):
        """Valid slices work like list slices."""
        batch = [{"id": 1}, {"id": 2}, {"id": 3}]
        results = fake_search(batch, total=3)
        for index in (
            slice(None),
            slice(None, None, -1),
            slice(1, None, -1),
            slice(5, 9),
            slice(-10, 10),
        ):
            with self.subTest(index=index):
                self.assertEqual(results[index], batch[index])
        empty = fake_search([], total=0)
        self.assertEqual(empty[:], [])
        self.assertEqual(empty[::-1], [])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
