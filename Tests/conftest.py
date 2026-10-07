# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""pytest configuration for the Biopython test suite.

The suite is run by pytest, configured here and in the
``[tool.pytest.ini_options]`` table of pyproject.toml.  ``run_tests.py`` is a
compatibility shim that translates its historical command line into a pytest
one.  This file keeps the behaviour of the bespoke runner that preceded it:

- ``--offline`` makes requires_internet skip and blocks non-local sockets.
- ``--check-skips`` fails the run if a module skips without being listed in
  Tests/expected_skips.txt.
- A module may skip only by raising MissingExternalDependencyError (or its
  subclass MissingPythonDependencyError) while being imported.  The same
  exceptions raised later, from code under test, are failures.
- Tests run from the Tests/ directory, with LANG restored before each module.
- A test module with no tests, or one that leaves the current directory
  changed, fails.
- The docstring examples of every Bio and BioSQL module run as doctests,
  collected by name through Tests/test_docstrings.py, one item per docstring.

Only unittest.TestCase subclasses are collected (python_classes and
python_functions are empty), exactly what run_tests.py used to run.
"""

import doctest
import gc
import importlib
import ipaddress
import os
import socket
import unittest
from fnmatch import fnmatchcase
from pathlib import Path
from pkgutil import iter_modules

import pytest

from Bio import MissingExternalDependencyError

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))

# These have been removed from Biopython and survive only as stubs whose
# import raises an informative ImportError pointing at the replacement.
EXCLUDE_DOCTEST_MODULES = [
    "Bio.Align.AlignInfo",
    "Bio.Align.Applications",
    "Bio.Alphabet",
    "Bio.Application",
    "Bio.Blast.Applications",
    "Bio.Emboss.Applications",
    "Bio.HMM",
    "Bio.pairwise2",
    "Bio.Phylo.Applications",
    "Bio.Sequencing.Applications",
]

# Exclude modules with online activity
# They are not excluded by default, use --offline to exclude them
ONLINE_DOCTEST_MODULES = [
    "Bio.Entrez",
    "Bio.ExPASy",
    "Bio.ExPASy.cellosaurus",
    "Bio.TogoWS",
    "Bio.UniProt",
]


try:
    import sqlite3

    del sqlite3
except ImportError:
    # May be missing on self-compiled Python
    EXCLUDE_DOCTEST_MODULES.append("Bio.SeqIO")
    EXCLUDE_DOCTEST_MODULES.append("Bio.SearchIO")


def find_modules(path):
    # Match setuptools.PackageFinder's built-in package exclusions.
    package_excludes = ("ez_setup", "*__pycache__")
    packages = set()
    for root, dirs, _ in os.walk(path, followlinks=True):
        candidates = dirs[:]
        dirs[:] = []
        for name in candidates:
            full_path = os.path.join(root, name)
            package = os.path.relpath(full_path, path).replace(os.path.sep, ".")
            # Skip directory trees that are not valid packages.
            if "." in name or not os.path.isfile(
                os.path.join(full_path, "__init__.py")
            ):
                continue
            if not any(fnmatchcase(package, pattern) for pattern in package_excludes):
                packages.add(package)
            # Keep searching subdirectories, as there may be more packages
            # down there, even if the parent was excluded.
            dirs.append(name)

    modules = set()
    for pkg in packages:
        modules.add(pkg)
        pkgpath = os.path.join(path, *pkg.split("."))
        for info in iter_modules([pkgpath]):
            if not info.ispkg:
                modules.add(pkg + "." + info.name)
    return modules


SYSTEM_LANG = os.environ.get("LANG", "C")  # Cache this

# Name of the test currently being run, used by the --offline socket guard
# to report which test attempted a network connection.
CURRENT_TEST = None


def _is_local_address(sock, address):
    """Check if a socket connection stays on this machine (PRIVATE).

    Allows AF_UNIX sockets and TCP/UDP connections to loopback addresses,
    e.g. for BioSQL tests run against a database server on this machine.
    """
    if sock.family == getattr(socket, "AF_UNIX", None):
        return True
    if isinstance(address, tuple) and address:
        host = address[0]
        if host == "localhost":
            return True
        try:
            return ipaddress.ip_address(host).is_loopback
        except ValueError:
            return False
    return False


def block_network_connections():
    """Patch socket.socket.connect to fail on non-local connections (PRIVATE).

    Monkeypatching urllib.request.urlopen is not enough: http.client,
    requests and raw sockets would still go online. Patching at the
    socket layer catches them all, whatever library made the attempt.
    """
    real_connect = socket.socket.connect

    def guarded_connect(sock, address):
        if _is_local_address(sock, address):
            return real_connect(sock, address)
        raise RuntimeError(
            f"{CURRENT_TEST or 'The test suite'} attempted a network "
            f"connection to {address!r} despite the --offline setting"
        )

    socket.socket.connect = guarded_connect


# Modules that skipped at import, mapped to the missing-dependency reason,
# keyed as in Tests/expected_skips.txt ("test_X" or "Bio.X").
_IMPORT_SKIPS = pytest.StashKey[dict]()
# Skips that --check-skips found missing from Tests/expected_skips.txt.
_UNEXPECTED_SKIPS = pytest.StashKey[list]()


def pytest_addoption(parser):
    group = parser.getgroup("biopython")
    group.addoption(
        "--offline",
        action="store_true",
        help="skip tests which require internet access, and fail any test "
        "that tries to make a non-local network connection",
    )
    group.addoption(
        "--check-skips",
        action="store_true",
        help="fail if a module skips that is not listed in Tests/expected_skips.txt",
    )


@pytest.hookimpl(trylast=True)
def pytest_configure(config):
    config.stash[_IMPORT_SKIPS] = {}
    if config.pluginmanager.has_plugin("dsession"):
        # This is the pytest-xdist controller, which runs no tests. Its
        # workers start in the invocation directory, resolve the command
        # line arguments there, and then change directory themselves.
        return
    if config.getoption("--offline"):
        # This is a bit of a hack...
        import requires_internet

        requires_internet.check.available = False
        # Block non-local socket connections so any test that tries
        # to use the internet fails loudly rather than going online.
        block_network_connections()
    # Always run tests from the Tests/ folder (as we assume this with
    # relative paths etc).  This runs after the built-in plugins have
    # resolved their own relative paths (--junitxml and so on) against the
    # invocation directory, and pytest resolves the test arguments against
    # that directory too, so pytest can be run from anywhere.
    os.chdir(TESTS_DIR)


def pytest_collectstart(collector):
    global CURRENT_TEST
    if isinstance(collector, pytest.Module):
        CURRENT_TEST = collector.nodeid


@pytest.hookimpl(tryfirst=True)
def pytest_runtest_setup(item):
    global CURRENT_TEST
    CURRENT_TEST = item.nodeid


def _skip_module(collector, name, error):
    """Skip a module whose import declared a missing dependency (PRIVATE)."""
    collector.config.stash[_IMPORT_SKIPS][name] = str(error)
    pytest.skip(str(error), allow_module_level=True)


class _TestModule(pytest.Module):
    """A Tests/test_*.py module (PRIVATE)."""

    def _getobj(self):
        try:
            return super()._getobj()
        except MissingExternalDependencyError as error:
            # Not an ImportError, so pytest lets it through unwrapped.
            _skip_module(self, self.path.stem, error)
        except self.CollectError as error:
            # pytest wraps an ImportError, which includes the subclass
            # MissingPythonDependencyError, in a CollectError.
            if isinstance(error.__cause__, MissingExternalDependencyError):
                _skip_module(self, self.path.stem, error.__cause__)
            raise

    def collect(self):
        tests = unittest.TestLoader().loadTestsFromModule(self.obj)
        if tests.countTestCases() == 0:
            raise self.CollectError(f"No tests found in {self.path.stem}")
        # Run the TestCase classes in unittest's order, sorted by name, as
        # run_tests.py did, not in pytest's definition order. Tests can
        # depend on it: under PyPy a file handle leaked by one class may only
        # be finalized during a later test's gc.collect().
        return sorted(super().collect(), key=lambda node: node.name)


class _DocTestRunner(doctest.DocTestRunner):
    """Hand doctest failures to pytest.DoctestItem (PRIVATE).

    DoctestItem.runtest() passes a list as ``out`` and reports every failure
    collected there with the expected and actual output, as pytest's own
    --doctest-modules does.  The remaining examples of a docstring still run
    after a failure, as they did under unittest's DocTestSuite.
    """

    def report_failure(self, out, test, example, got):
        out.append(doctest.DocTestFailure(test, example, got))

    def report_unexpected_exception(self, out, test, example, exc_info):
        out.append(doctest.UnexpectedException(test, example, exc_info))


class _DoctestItem(pytest.DoctestItem):
    """The examples in one docstring (PRIVATE)."""

    def reportinfo(self):
        # Point at the docstring in the Bio source, not Tests/test_docstrings.py.
        return Path(self.dtest.filename), self.dtest.lineno, f"[doctest] {self.name}"


class _DocstringModule(pytest.Module):
    """The docstring examples of one Bio or BioSQL module (PRIVATE).

    The module is imported by name rather than from a path, so the tests use
    whichever Bio is on sys.path: the source tree, or an installed copy.
    """

    def _getobj(self):
        try:
            return importlib.import_module(self.name)
        except MissingExternalDependencyError as error:
            _skip_module(self, self.name, error)

    def collect(self):
        module = self.obj
        runner = _DocTestRunner(verbose=False, optionflags=doctest.ELLIPSIS)
        for test in sorted(doctest.DocTestFinder().find(module)):
            if not test.examples:
                continue
            if not test.filename:
                test.filename = module.__file__
            yield _DoctestItem.from_parent(
                self, name=test.name, runner=runner, dtest=test
            )


class _DocstringModules(pytest.File):
    """Collect the doctests of every Bio and BioSQL module (PRIVATE).

    This is anchored at Tests/test_docstrings.py, so it runs with the rest of
    the suite and ``pytest test_docstrings.py::Bio.Seq`` selects one module.
    """

    def collect(self):
        excluded = set(EXCLUDE_DOCTEST_MODULES)
        if self.config.getoption("--offline"):
            excluded.update(ONLINE_DOCTEST_MODULES)
        names = find_modules(os.path.join(TESTS_DIR, os.pardir)) - excluded
        for name in sorted(names):
            yield _DocstringModule.from_parent(
                self, path=self.path, name=name, nodeid=f"{self.nodeid}::{name}"
            )


def pytest_pycollect_makemodule(module_path, parent):
    if module_path.name == "test_docstrings.py":
        return _DocstringModules.from_parent(parent, path=module_path)
    return _TestModule.from_parent(parent, path=module_path)


@pytest.fixture(autouse=True, scope="module")
def _module_hygiene():
    """Give each module the clean start run_tests.py used to give it (PRIVATE)."""
    # Restore the language and thus default encoding (in case a prior
    # test changed this, e.g. to help with detecting command line tools)
    os.environ["LANG"] = SYSTEM_LANG
    cwd = os.getcwd()
    yield
    # Running under PyPy we were leaking file handles...
    gc.collect()
    now = os.getcwd()
    if now != cwd:
        os.chdir(cwd)
        pytest.fail(f"Current directory changed\nWas: {cwd}\nNow: {now}", pytrace=False)


def pytest_sessionfinish(session):
    config = session.config
    if not config.getoption("--check-skips"):
        return
    expected = set()
    with open(os.path.join(TESTS_DIR, "expected_skips.txt")) as handle:
        for line in handle:
            entry = line.split("#", 1)[0].strip()
            if entry:
                expected.add(entry)
    unexpected = sorted(set(config.stash[_IMPORT_SKIPS]) - expected)
    config.stash[_UNEXPECTED_SKIPS] = unexpected
    if unexpected and session.exitstatus == pytest.ExitCode.OK:
        session.exitstatus = pytest.ExitCode.TESTS_FAILED


def pytest_terminal_summary(terminalreporter, config):
    skips = config.stash.get(_IMPORT_SKIPS, {})
    if skips:
        terminalreporter.section("modules skipped at import")
        for name in sorted(skips):
            terminalreporter.line(f"{name} -- {skips[name]}")
    unexpected = config.stash.get(_UNEXPECTED_SKIPS, [])
    if unexpected:
        terminalreporter.section("unexpected skips", red=True)
        terminalreporter.line(
            f"FAILED (unexpected skips = {len(unexpected)})\n"
            "These modules skipped but are not listed in "
            "Tests/expected_skips.txt:",
            red=True,
        )
        for name in unexpected:
            terminalreporter.line(f"    {name} -- {skips[name]}")
