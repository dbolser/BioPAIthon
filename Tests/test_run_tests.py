# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for the Biopython test runner."""

import os
import subprocess
import sys
import tempfile
import textwrap
import unittest

from conftest import find_modules
from run_tests import pytest_args

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
PYPROJECT = os.path.join(TESTS_DIR, os.pardir, "pyproject.toml")


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
                "ez_setup/sub/__init__.py",
                "ez_setup/sub/mod.py",
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
                    # An excluded package is not itself reported, but it is
                    # still descended into, exactly as setuptools does.
                    "ez_setup.sub",
                    "ez_setup.sub.mod",
                },
            )


class ShimTests(unittest.TestCase):
    """run_tests.py translates its historical arguments for pytest."""

    def path(self, name):
        return os.path.join(TESTS_DIR, name)

    def test_no_arguments_runs_everything(self):
        self.assertEqual(pytest_args([]), [TESTS_DIR])
        self.assertEqual(pytest_args(["--offline"]), ["--offline", TESTS_DIR])

    def test_options(self):
        self.assertEqual(
            pytest_args(["--offline", "--check-skips", "-v", "--verbose", "-x"]),
            ["--offline", "--check-skips", "-v", "--verbose", "-x", TESTS_DIR],
        )

    def test_test_modules(self):
        self.assertEqual(
            pytest_args(["--offline", "test_Seq_objs", "test_GenBank.py"]),
            [
                "--offline",
                self.path("test_Seq_objs.py"),
                self.path("test_GenBank.py"),
            ],
        )

    def test_doctests(self):
        docstrings = self.path("test_docstrings.py")
        self.assertEqual(pytest_args(["doctest"]), [docstrings])
        self.assertEqual(
            pytest_args(["Bio.Seq", "BioSQL.BioSeq"]),
            [docstrings + "::Bio.Seq", docstrings + "::BioSQL.BioSeq"],
        )

    def test_legacy_options_are_ignored(self):
        self.assertEqual(
            pytest_args(["--doctest", "-g", "--generate", "test_Seq_objs"]),
            [self.path("test_Seq_objs.py")],
        )

    def test_option_values_pass_through(self):
        # A value that looks like a test name is left alone after -k, and
        # the value of a long option does not count as a test selection.
        self.assertEqual(
            pytest_args(["-k", "test_translation", "--tb", "short", "-n", "4"]),
            ["-k", "test_translation", "--tb", "short", "-n", "4", TESTS_DIR],
        )

    def test_node_ids_pass_through(self):
        node = "test_Seq_objs.py::StringMethodTests"
        self.assertEqual(pytest_args([node]), [node])

    def test_help_without_site_packages(self):
        """--help must not need pytest, setuptools or any other site package."""
        result = subprocess.run(
            [sys.executable, "-S", "run_tests.py", "--help"],
            cwd=TESTS_DIR,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--offline", result.stdout)

    def test_missing_pytest_explains(self):
        """Without pytest the shim exits 2 with an install hint."""
        probe = subprocess.run(
            [sys.executable, "-S", "-c", "import pytest"], capture_output=True
        )
        if probe.returncode == 0:
            self.skipTest("pytest is importable without site-packages")
        result = subprocess.run(
            [sys.executable, "-S", "run_tests.py", "--offline", "test_Seq_objs"],
            cwd=TESTS_DIR,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertIn("pip install pytest", result.stderr)


class ProbeTests(unittest.TestCase):
    """Run pytest in a subprocess, with Tests/conftest.py, on throwaway modules.

    These pin the behaviour Tests/conftest.py builds on pytest internals, so
    that a pytest release which changes them fails here rather than quietly
    changing what the suite does.
    """

    def run_pytest(self, files, *options):
        """Write files to a temporary directory and run pytest on them."""
        with tempfile.TemporaryDirectory() as directory:
            paths = []
            for name, text in files.items():
                path = os.path.join(directory, name)
                with open(path, "w") as handle:
                    handle.write(textwrap.dedent(text))
                if name.startswith("test_"):
                    paths.append(path)
            env = dict(os.environ)
            env["PYTHONPATH"] = os.pathsep.join(
                [directory, *filter(None, [env.get("PYTHONPATH")])]
            )
            # --noconftest and -c keep any conftest.py or configuration file
            # in the temporary directory's parents out of the run; -p loads
            # Tests/conftest.py (importable from the cwd) as a plugin.
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "pytest",
                    "--noconftest",
                    "-p",
                    "conftest",
                    "-p",
                    "no:cacheprovider",
                    "-c",
                    PYPROJECT,
                    "--rootdir",
                    directory,
                    *options,
                    *paths,
                ],
                cwd=TESTS_DIR,
                env=env,
                capture_output=True,
                text=True,
            )
        return result.returncode, result.stdout + result.stderr

    def test_import_time_skips_and_check_skips(self):
        returncode, output = self.run_pytest(
            {
                "test_probe_external.py": """\
                    from Bio import MissingExternalDependencyError
                    raise MissingExternalDependencyError("no probe tool")
                """,
                "test_probe_python.py": """\
                    from Bio import MissingPythonDependencyError
                    raise MissingPythonDependencyError("no probe package")
                """,
                "test_probe_ok.py": """\
                    import unittest
                    class Ok(unittest.TestCase):
                        def test_ok(self):
                            pass
                """,
            },
            "--check-skips",
        )
        self.assertEqual(returncode, 1, output)
        self.assertIn("1 passed, 2 skipped", output)
        self.assertIn("FAILED (unexpected skips = 2)", output)
        self.assertIn("test_probe_external -- no probe tool", output)
        self.assertIn("test_probe_python -- no probe package", output)

    def test_listed_skips_pass_check_skips(self):
        # Borrow the name of a test module listed in Tests/expected_skips.txt.
        with open(os.path.join(TESTS_DIR, "expected_skips.txt")) as handle:
            listed = next(line.split()[0] for line in handle if line[:5] == "test_")
        returncode, output = self.run_pytest(
            {
                f"{listed}.py": """\
                    from Bio import MissingExternalDependencyError
                    raise MissingExternalDependencyError("no probe tool")
                """,
                "test_probe_ok.py": """\
                    import unittest
                    class Ok(unittest.TestCase):
                        def test_ok(self):
                            pass
                """,
            },
            "--check-skips",
        )
        self.assertEqual(returncode, 0, output)
        self.assertIn("1 passed, 1 skipped", output)
        self.assertNotIn("unexpected skips", output)

    def test_failures(self):
        """Late dependency errors, import errors, empty modules and chdir fail."""
        returncode, output = self.run_pytest(
            {
                "test_probe_late.py": """\
                    import unittest
                    from Bio import MissingExternalDependencyError
                    class Late(unittest.TestCase):
                        def test_late(self):
                            raise MissingExternalDependencyError("too late")
                """,
                "test_probe_import.py": """\
                    import probe_module_that_does_not_exist
                """,
                "test_probe_empty.py": """\
                    import unittest
                    class Empty(unittest.TestCase):
                        pass
                """,
                "test_probe_chdir.py": """\
                    import os
                    import unittest
                    class Chdir(unittest.TestCase):
                        def test_chdir(self):
                            os.chdir(os.pardir)
                """,
                "test_probe_ok.py": """\
                    import unittest
                    class Ok(unittest.TestCase):
                        def test_ok(self):
                            pass
                """,
            }
        )
        self.assertEqual(returncode, 1, output)
        # The collection errors do not stop the other modules running.
        self.assertIn("1 failed, 2 passed, 3 errors", output)
        self.assertRegex(output, r"FAILED \S*test_probe_late.py::Late::test_late")
        self.assertRegex(output, r"ERROR \S*test_probe_import.py\n")
        self.assertIn("No tests found in test_probe_empty", output)
        self.assertRegex(output, r"ERROR \S*test_probe_chdir.py::Chdir::test_chdir")
        self.assertIn("Current directory changed", output)

    def test_classes_run_in_unittest_order(self):
        returncode, output = self.run_pytest(
            {
                "test_probe_order.py": """\
                    import unittest
                    class Second(unittest.TestCase):
                        def test_b(self):
                            pass
                        def test_a(self):
                            pass
                    class First(unittest.TestCase):
                        def test_z(self):
                            pass
                """,
            },
            "--collect-only",
            "-q",
        )
        self.assertEqual(returncode, 0, output)
        self.assertRegex(
            output,
            r"test_probe_order.py::First::test_z\n"
            r"\S*test_probe_order.py::Second::test_a\n"
            r"\S*test_probe_order.py::Second::test_b\n",
        )

    def test_offline_blocks_network(self):
        returncode, output = self.run_pytest(
            {
                "test_probe_network.py": """\
                    import socket
                    import unittest
                    class Network(unittest.TestCase):
                        def test_connect(self):
                            # TEST-NET-1, so no DNS lookup is needed.
                            socket.create_connection(("192.0.2.1", 80), timeout=5)
                """,
            },
            "--offline",
        )
        self.assertEqual(returncode, 1, output)
        self.assertIn(
            "test_probe_network.py::Network::test_connect attempted a network "
            "connection to ('192.0.2.1', 80) despite the --offline setting",
            output,
        )

    def test_doctest_failure_is_reported_inline(self):
        # The probe plugin anchors the doctest collector at a test file, as
        # Tests/test_docstrings.py does, but for a throwaway module.
        returncode, output = self.run_pytest(
            {
                "probe_docstrings.py": '''\
                    """A module with one broken docstring example."""


                    def add():
                        """Add one and one.

                        >>> 1 + 1
                        3
                        """
                ''',
                "probe_plugin.py": """\
                    import pytest
                    from conftest import _DocstringModule

                    class Anchor(pytest.File):
                        def collect(self):
                            yield _DocstringModule.from_parent(
                                self,
                                path=self.path,
                                name="probe_docstrings",
                                nodeid=self.nodeid + "::probe_docstrings",
                            )

                    @pytest.hookimpl(tryfirst=True)
                    def pytest_pycollect_makemodule(module_path, parent):
                        if module_path.name == "test_probe_anchor.py":
                            return Anchor.from_parent(parent, path=module_path)
                """,
                "test_probe_anchor.py": "",
            },
            "-p",
            "probe_plugin",
            "-vv",
        )
        self.assertEqual(returncode, 1, output)
        self.assertIn("1 failed", output)
        # Verbose output and the failure both point at the docstring's source.
        self.assertRegex(
            output,
            r"test_probe_anchor.py::probe_docstrings::probe_docstrings.add"
            r" <- \S*probe_docstrings.py FAILED",
        )
        self.assertRegex(output, r"probe_docstrings.py:7: DocTestFailure")
        # Python 3.13 and later strip the docstring's indentation.
        self.assertRegex(output, r"007 +>>> 1 \+ 1")
        self.assertIn("Expected:\n    3\nGot:\n    2", output)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
