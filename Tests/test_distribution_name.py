"""Tests for answering importlib.metadata queries for "biopython".

The distribution is called biopaithon, the package is still Bio, and code in
the wild asks importlib.metadata for the version of "biopython".
"""

import importlib.metadata
import sys
import unittest

import Bio


def _installed(name):
    try:
        importlib.metadata.distribution(name)
    except importlib.metadata.PackageNotFoundError:
        return False
    return True


class TestDistributionName(unittest.TestCase):
    def setUp(self):
        finders = [
            f for f in sys.meta_path if isinstance(f, Bio._BiopythonDistributionFinder)
        ]
        self.assertEqual(len(finders), 1, "exactly one finder is registered")
        self.finder = finders[0]

    def test_registration_is_idempotent(self):
        Bio._register_distribution_finder()
        finders = [
            f for f in sys.meta_path if isinstance(f, Bio._BiopythonDistributionFinder)
        ]
        self.assertEqual(len(finders), 1)

    def test_other_names_are_not_answered(self):
        Context = importlib.metadata.DistributionFinder.Context
        for name in (None, "numpy", "biopaithon", "bio"):
            with self.subTest(name=name):
                self.assertEqual(
                    list(self.finder.find_distributions(Context(name=name))), []
                )

    @unittest.skipUnless(_installed("biopaithon"), "biopaithon is not installed")
    def test_biopython_resolves_to_biopaithon(self):
        expected = importlib.metadata.version("biopaithon")
        for name in ("biopython", "Biopython", "BIOPYTHON"):
            with self.subTest(name=name):
                self.assertEqual(importlib.metadata.version(name), expected)
        self.assertEqual(
            importlib.metadata.distribution("biopython").metadata["Name"], "biopaithon"
        )

    @unittest.skipUnless(_installed("biopaithon"), "biopaithon is not installed")
    def test_not_listed_twice(self):
        names = [d.metadata["Name"] for d in importlib.metadata.distributions()]
        self.assertEqual(names.count("biopaithon"), 1)

    @unittest.skipIf(_installed("biopaithon"), "biopaithon is installed")
    def test_nothing_to_answer_with(self):
        with self.assertRaises(importlib.metadata.PackageNotFoundError):
            importlib.metadata.version("biopython")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
