"""Tests for answering importlib.metadata queries for "biopython".

The distribution is called biopaithon, the package is still Bio, and code in
the wild asks importlib.metadata for the version of "biopython".

The tests that resolve names see only metadata they fabricate themselves, so
they do not depend on what the environment running them holds: an installed
biopaithon or none, a real Biopython or none, or the stale biopython.egg-info a
checkout keeps from before the rename.
"""

import contextlib
import importlib.metadata
import json
import os
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import Bio

try:
    import tomllib
except ImportError:  # Python 3.10
    tomllib = None

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
PYPROJECT = os.path.join(ROOT, "pyproject.toml")

# Versions no real release has, so that an answer from the environment rather
# than from the fabricated metadata cannot pass unnoticed.
BIOPAITHON = ("biopaithon", "1.89.dev0+fixture")
BIOPYTHON = ("biopython", "1.85+fixture")


class TestDistributionName(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # sys.path as a bare interpreter has it: the standard library, which
        # importlib.metadata imports from as it goes, and no site-packages or
        # working directory that could hold another biopython or biopaithon.
        code = "import json, sys; print(json.dumps(sys.path))"
        output = subprocess.check_output([sys.executable, "-I", "-S", "-c", code])
        cls.stdlib_path = json.loads(output)

    def setUp(self):
        finders = [
            f for f in sys.meta_path if isinstance(f, Bio._BiopythonDistributionFinder)
        ]
        self.assertEqual(len(finders), 1, "exactly one finder is registered")
        self.finder = finders[0]

    @contextlib.contextmanager
    def installed_only(self, *distributions):
        """Make importlib.metadata see only these (name, version) distributions.

        Each gets the dist-info directory an installer would write, in a
        temporary directory that replaces everything on sys.path meanwhile but
        the standard library.
        """
        with tempfile.TemporaryDirectory() as directory:
            for name, version in distributions:
                dist_info = os.path.join(directory, f"{name}-{version}.dist-info")
                os.mkdir(dist_info)
                with open(os.path.join(dist_info, "METADATA"), "w") as handle:
                    handle.write(
                        f"Metadata-Version: 2.1\nName: {name}\nVersion: {version}\n"
                    )
            with mock.patch.object(sys, "path", [directory, *self.stdlib_path]):
                yield

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

    def test_biopython_resolves_to_biopaithon(self):
        names = ("biopython", "Biopython", "BIOPYTHON")
        with self.installed_only(BIOPAITHON):
            versions = {name: importlib.metadata.version(name) for name in names}
            answer = importlib.metadata.distribution("biopython").metadata["Name"]
        self.assertEqual(versions, dict.fromkeys(names, BIOPAITHON[1]))
        self.assertEqual(answer, "biopaithon")

    def test_real_biopython_answers_first(self):
        with self.installed_only(BIOPYTHON, BIOPAITHON):
            version = importlib.metadata.version("biopython")
            answer = importlib.metadata.distribution("biopython").metadata["Name"]
        self.assertEqual(version, BIOPYTHON[1])
        self.assertEqual(answer, "biopython")

    def test_nothing_to_answer_with(self):
        with self.installed_only():
            with self.assertRaises(importlib.metadata.PackageNotFoundError):
                importlib.metadata.version("biopython")
            # A search finds nothing, rather than failing from inside the finder.
            found = list(importlib.metadata.distributions(name="biopython"))
        self.assertEqual(found, [])

    def test_not_listed_twice(self):
        for installed in [(BIOPAITHON,), (BIOPYTHON, BIOPAITHON)]:
            with self.subTest(installed=installed):
                with self.installed_only(*installed):
                    dists = importlib.metadata.distributions()
                    names = [dist.metadata["Name"] for dist in dists]
                # The standard library may hold records of its own (PyPy's
                # does), so count only the two names in question.
                self.assertEqual(names.count("biopaithon"), 1)
                self.assertEqual(names.count("biopython"), installed.count(BIOPYTHON))

    @unittest.skipIf(tomllib is None, "needs tomllib, new in Python 3.11")
    def test_finder_names_this_project(self):
        # The tests above fabricate "biopaithon", the name the finder looks up,
        # so they would still pass if the project were renamed and the finder
        # were not.
        with open(PYPROJECT, "rb") as handle:
            name = tomllib.load(handle)["project"]["name"]
        self.assertEqual(name, BIOPAITHON[0])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
