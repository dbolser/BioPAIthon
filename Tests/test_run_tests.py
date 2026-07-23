# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for the Biopython test runner."""

import os
import subprocess
import sys
import tempfile
import unittest

from run_tests import find_modules


class FindModulesTests(unittest.TestCase):
    def test_classic_packages_and_modules(self):
        with tempfile.TemporaryDirectory() as directory:
            files = [
                "example/__init__.py",
                "example/module.py",
                "example/subpackage/__init__.py",
                "example/subpackage/child.py",
                "example/ez_setup/__init__.py",
                "example/ez_setup/helper.py",
                "example/data/not_a_module.py",
                "example/dotted.name/__init__.py",
                "example/__pycache__/__init__.py",
                "example/project__pycache__/__init__.py",
                "not_a_package/nested/__init__.py",
                "ez_setup/__init__.py",
            ]
            for filename in files:
                path = os.path.join(directory, filename)
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with open(path, "w"):
                    pass

            self.assertEqual(
                find_modules(directory),
                {
                    "example",
                    "example.ez_setup",
                    "example.ez_setup.helper",
                    "example.module",
                    "example.subpackage",
                    "example.subpackage.child",
                },
            )

    def test_runner_starts_without_site_packages(self):
        """The test runner must not require setuptools or another site package."""
        result = subprocess.run(
            [sys.executable, "-S", "run_tests.py", "--help"],
            cwd=os.path.dirname(__file__),
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--offline", result.stdout)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
