#!/usr/bin/env python
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Run a set of PyUnit-based regression tests.

This will find all modules whose name is "test_*.py" in the test
directory, and run them.  Various command line options provide
additional facilities.

Command line options::

    --help        -- show usage info
    --offline     -- skip tests which require internet access
    --check-skips -- fail if a module skips that is not listed in
                     Tests/expected_skips.txt
    -v;--verbose  -- run tests with higher verbosity
    <test_name>   -- supply the name of one (or more) tests to be run.
                     The .py file extension is optional.
    doctest       -- run the docstring tests.

<test_name> and doctest, if supplied, must come after the other options.
By default, all tests are run.
"""

# standard modules
import doctest
import gc
import getopt
import ipaddress
import os
import socket
import sys
import time
import traceback
import unittest
from fnmatch import fnmatchcase
from io import StringIO
from pkgutil import iter_modules

# The default verbosity (not verbose)
VERBOSITY = 0

# These have been removed from Biopython and survive only as stubs whose
# import raises an informative ImportError pointing at the replacement.
EXCLUDE_DOCTEST_MODULES = [
    "Bio.Align.Applications",
    "Bio.Alphabet",
    "Bio.Application",
    "Bio.Blast.Applications",
    "Bio.Emboss.Applications",
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


def main(argv):
    """Run tests, return number of failures (integer)."""
    # Using "export LANG=C" (which should work on Linux and similar) can
    # avoid problems detecting optional command line tools on
    # non-English OS (we may want 'command not found' in English).
    # HOWEVER, we do not want to change the default encoding which is
    # rather important on Python 3 with unicode.
    # lang = os.environ['LANG']

    # get the command line options
    try:
        opts, args = getopt.getopt(
            argv,
            "gv",
            ["generate", "verbose", "doctest", "help", "offline", "check-skips"],
        )
    except getopt.error as msg:
        print(msg)
        print(__doc__)
        return 2

    verbosity = VERBOSITY
    check_skips = False

    # deal with the options
    for opt, _ in opts:
        if opt == "--help":
            print(__doc__)
            return 0
        if opt == "--check-skips":
            check_skips = True
        if opt == "--offline":
            print("Skipping any tests requiring internet access")
            EXCLUDE_DOCTEST_MODULES.extend(ONLINE_DOCTEST_MODULES)
            # This is a bit of a hack...
            import requires_internet

            requires_internet.check.available = False
            # Block non-local socket connections so any test that tries
            # to use the internet fails loudly rather than going online.
            block_network_connections()

        if opt == "-v" or opt == "--verbose":
            verbosity = 2

    # deal with the arguments, which should be names of tests to run
    for arg_num in range(len(args)):
        # strip off the .py if it was included
        if args[arg_num][-3:] == ".py":
            args[arg_num] = args[arg_num][:-3]

    print(f"Python version: {sys.version}")
    print(f"Operating system: {os.name} {sys.platform}")

    # run the tests
    runner = TestRunner(args, verbosity, check_skips)
    return runner.run()


class TestRunner(unittest.TextTestRunner):
    if __name__ == "__main__":
        file = sys.argv[0]
    else:
        file = __file__
    testdir = os.path.abspath(os.path.dirname(file) or os.curdir)

    def __init__(self, tests=None, verbosity=0, check_skips=False):
        """Initialise test runner.

        If not tests are specified, we run them all,
        including the doctests.

        Defaults to running without any verbose logging.
        """
        if tests is None:
            self.tests = []
        else:
            self.tests = tests
        self.check_skips = check_skips
        # Modules which skipped, mapped to the missing-dependency reason.
        self.skips = {}
        # Number of individual test cases run across all modules.
        self.cases_run = 0
        if not self.tests:
            # Make a list of all applicable test modules.
            names = os.listdir(TestRunner.testdir)
            for name in names:
                if name[:5] == "test_" and name[-3:] == ".py":
                    self.tests.append(name[:-3])
            self.tests.sort()
            self.tests.append("doctest")
        if "doctest" in self.tests:
            self.tests.remove("doctest")
            modules = find_modules(self.testdir + "/..")
            modules.difference_update(set(EXCLUDE_DOCTEST_MODULES))
            self.tests.extend(sorted(modules))
        stream = StringIO()
        unittest.TextTestRunner.__init__(self, stream, verbosity=verbosity)

    def runTest(self, name):
        """Run one test module; return "ok", "skip" or "fail" (PRIVATE).

        A module may skip only by raising MissingExternalDependencyError
        (or its subclass MissingPythonDependencyError) while being
        imported -- that is a test module declaring up front that a
        dependency is absent.  The same exceptions raised later, from
        code under test, are real failures and are reported as such.
        """
        # MissingPythonDependencyError subclasses this; one except covers both.
        from Bio import MissingExternalDependencyError

        global CURRENT_TEST
        CURRENT_TEST = name
        result = self._makeResult()
        output = StringIO()
        # Restore the language and thus default encoding (in case a prior
        # test changed this, e.g. to help with detecting command line tools)
        os.environ["LANG"] = SYSTEM_LANG
        # Always run tests from the Tests/ folder where run_tests.py
        # should be located (as we assume this with relative paths etc)
        os.chdir(self.testdir)
        try:
            stdout = sys.stdout
            sys.stdout = output
            if name.startswith("test_"):
                # It's a unittest
                sys.stderr.write(f"{name} ... ")
                try:
                    module = __import__(name)
                except MissingExternalDependencyError as msg:
                    self.skips[name] = str(msg)
                    sys.stderr.write(f"skipping. {msg}\n")
                    return "skip"
                loader = unittest.TestLoader()
                suite = loader.loadTestsFromModule(module)
                if loader.errors:
                    sys.stderr.write("loading tests failed:\n")
                    for msg in loader.errors:
                        sys.stderr.write(f"{msg}\n")
                    return "fail"
                if suite.countTestCases() == 0:
                    raise RuntimeError(f"No tests found in {name}")
            else:
                # It's a doc test
                sys.stderr.write(f"{name} docstring test ... ")
                try:
                    module = __import__(name, fromlist=name.split("."))
                except MissingExternalDependencyError as msg:
                    self.skips[name] = str(msg)
                    sys.stderr.write(f"skipping. {msg}\n")
                    return "skip"
                except ImportError as e:
                    sys.stderr.write("FAIL, ImportError\n")
                    result.stream.write(f"ERROR while importing {name}: {e}\n")
                    result.printErrors()
                    return "fail"
                suite = doctest.DocTestSuite(module, optionflags=doctest.ELLIPSIS)
                del module
            suite.run(result)
            self.cases_run += result.testsRun
            if self.testdir != os.path.abspath("."):
                sys.stderr.write("FAIL\n")
                result.stream.write(result.separator1 + "\n")
                result.stream.write(f"ERROR: {name}\n")
                result.stream.write(result.separator2 + "\n")
                result.stream.write("Current directory changed\n")
                result.stream.write(f"Was: {self.testdir}\n")
                result.stream.write(f"Now: {os.path.abspath('.')}\n")
                os.chdir(self.testdir)
                if not result.wasSuccessful():
                    result.printErrors()
                return "fail"
            elif result.wasSuccessful():
                sys.stderr.write("ok\n")
                return "ok"
            else:
                sys.stderr.write("FAIL\n")
                result.printErrors()
            return "fail"
        except Exception:  # noqa: BLE001
            # An unexpected error, e.g. during import or test loading.
            # MissingExternalDependencyError raised from code under test
            # (rather than while importing the test module) also lands
            # here, deliberately: that is a failure, not a skip.
            sys.stderr.write("ERROR\n")
            result.stream.write(result.separator1 + "\n")
            result.stream.write(f"ERROR: {name}\n")
            result.stream.write(result.separator2 + "\n")
            result.stream.write(traceback.format_exc())
            return "fail"
        finally:
            sys.stdout = stdout
            # Running under PyPy we were leaking file handles...
            gc.collect()

    def checkSkips(self):
        """Report skips absent from the manifest; return their number (PRIVATE)."""
        manifest = os.path.join(self.testdir, "expected_skips.txt")
        expected = set()
        with open(manifest) as handle:
            for line in handle:
                entry = line.split("#", 1)[0].strip()
                if entry:
                    expected.add(entry)
        unexpected = sorted(set(self.skips) - expected)
        if unexpected:
            sys.stderr.write(
                "FAILED (unexpected skips = %d)\n"
                "These modules skipped but are not listed in "
                "Tests/expected_skips.txt:\n" % len(unexpected)
            )
            for name in unexpected:
                sys.stderr.write(f"    {name} -- {self.skips[name]}\n")
        return len(unexpected)

    def run(self):
        """Run tests, return number of failures (integer)."""
        failures = 0
        start_time = time.time()
        for test in self.tests:
            status = self.runTest(test)
            if status == "fail":
                failures += 1
        total = len(self.tests)
        skips = len(self.skips)
        stop_time = time.time()
        time_taken = stop_time - start_time
        sys.stderr.write(self.stream.getvalue())
        sys.stderr.write("-" * 70 + "\n")
        sys.stderr.write(
            "Ran %d module%s (%d case%s) in %.3f seconds, %d skipped, %d failed\n"
            % (
                total,
                "s" if total != 1 else "",
                self.cases_run,
                "s" if self.cases_run != 1 else "",
                time_taken,
                skips,
                failures,
            )
        )
        sys.stderr.write("\n")
        if failures:
            sys.stderr.write("FAILED (failures = %d)\n" % failures)
        if self.check_skips:
            failures += self.checkSkips()
        return failures


if __name__ == "__main__":
    errors = main(sys.argv[1:])
    if errors:
        # Doing a sys.exit(...) isn't nice if run from IDLE...
        sys.exit(1)
