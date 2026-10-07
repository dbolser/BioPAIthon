#!/usr/bin/env python
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Run the Biopython test suite (a compatibility shim for pytest).

The suite is run by pytest, configured in pyproject.toml and
Tests/conftest.py.  This script translates its historical command line into
a pytest one, so ``python run_tests.py --offline`` keeps working.

Command line options::

    --help        -- show usage info
    --offline     -- skip tests which require internet access
    --check-skips -- fail if a module skips that is not listed in
                     Tests/expected_skips.txt
    -v;--verbose  -- run tests with higher verbosity
    <test_name>   -- supply the name of one (or more) tests to be run.
                     The .py file extension is optional.
    doctest       -- run the docstring tests.
    Bio.X         -- run the docstring tests of one module.

By default, all tests are run, including the docstring tests.  Any other
argument is passed to pytest unchanged, e.g. -x, -k EXPRESSION or a node ID
such as test_Seq_objs.py::StringMethodTests.  The legacy --doctest and
-g/--generate options are accepted and ignored.

The pytest equivalent of ``python run_tests.py --offline`` is
``python -m pytest --offline``, which also works from the repository root.
"""

import os
import re
import sys

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))

# pytest's single-letter options that take a value as the next argument,
# which must not be mistaken for a test name (as in -k test_translation).
VALUE_OPTIONS = {"-c", "-k", "-m", "-n", "-o", "-p", "-r", "-W"}

IGNORED_OPTIONS = {"--doctest", "-g", "--generate"}

INSTALL_HINT = """\
run_tests.py now runs the test suite with pytest, which is not installed.
Install it with:  python -m pip install pytest
or install all the test dependencies with:  python -m pip install -e ".[test]"
"""


def pytest_args(argv):
    """Translate run_tests.py arguments into pytest arguments (PRIVATE)."""
    args = []
    selected = False
    expect_value = False
    for arg in argv:
        if expect_value:
            args.append(arg)
            expect_value = False
        elif arg in IGNORED_OPTIONS:
            print(f"Ignoring {arg}, which no longer has any effect.")
        elif arg.startswith("-"):
            # Including --offline, --check-skips and -v/--verbose, which
            # mean the same to pytest (see Tests/conftest.py).
            args.append(arg)
            expect_value = arg in VALUE_OPTIONS
        elif arg == "doctest":
            args.append(os.path.join(TESTS_DIR, "test_docstrings.py"))
            selected = True
        elif re.fullmatch(r"test_\w+(\.py)?", arg):
            name = arg.removesuffix(".py")
            args.append(os.path.join(TESTS_DIR, name + ".py"))
            selected = True
        elif re.fullmatch(r"(Bio|BioSQL)(\.\w+)*", arg):
            args.append(os.path.join(TESTS_DIR, "test_docstrings.py") + "::" + arg)
            selected = True
        else:
            # Anything else goes to pytest unchanged: a node ID, a path, or
            # the value of a long option, as in --tb short.
            args.append(arg)
            selected = selected or "::" in arg or os.path.exists(arg)
    if not selected:
        args.append(TESTS_DIR)
    return args


def main(argv):
    """Run the tests and return pytest's exit code."""
    if "--help" in argv:
        print(__doc__)
        return 0
    try:
        import pytest
    except ImportError:
        sys.stderr.write(INSTALL_HINT)
        return 2
    return pytest.main(pytest_args(argv))


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
