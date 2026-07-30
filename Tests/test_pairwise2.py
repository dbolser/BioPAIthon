# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for pairwise2 module using the default C functions.

This test file imports the TestCases from ``pairwise2_testCases.py``.
If you want to add more tests, do this over there.

"""

import platform
import subprocess
import sys
import sysconfig
import unittest

# Import all test classes deliberately:
from pairwise2_testCases import *  # noqa: F401, F403

# Implicitly using functions from C extension:
from Bio import pairwise2

if pairwise2.rint == pairwise2._python_rint:
    from Bio import MissingExternalDependencyError

    raise MissingExternalDependencyError("Missing or non-compiled file: 'cpairwise2'")


class RintTest(unittest.TestCase):
    def test_precision(self):
        self.assertEqual(pairwise2.rint(1.25), 1250)
        self.assertEqual(pairwise2.rint(1.25, 100), 125)
        self.assertEqual(pairwise2.rint(1.25, precision=100), 125)

    def test_precision_int_boundaries(self):
        self.assertEqual(pairwise2.rint(0.0, 2**31 - 1), 0)
        self.assertEqual(pairwise2.rint(0.0, -(2**31)), 0)

    def test_precision_overflow(self):
        for precision in (2**31, -(2**31) - 1, 2**62):
            with self.subTest(precision=precision):
                with self.assertRaises(OverflowError):
                    pairwise2.rint(0.0, precision)


_stable_cpython_refcounts = (
    platform.python_implementation() == "CPython"
    and not sysconfig.get_config_var("Py_GIL_DISABLED")
)

# One sequence that converts to ASCII bytes and one that cannot. This is the
# only combination that reaches the conversion failure path in
# _make_score_matrix_fast while leaving a live object behind to mishandle;
# when neither sequence converts there is nothing to release twice, which is
# why the crash looked intermittent to the reporter of biopython/biopython#3771.
_CONVERTIBLE = b"ACGTACGTACGT"
_UNCONVERTIBLE = "üüü"


class BytesConversionOwnershipTest(unittest.TestCase):
    """Reference handling on the sequence-to-bytes conversion path."""

    @unittest.skipUnless(
        _stable_cpython_refcounts,
        "GIL-enabled CPython reference counts are required",
    )
    def test_bytes_sequence_reference_not_stolen(self):
        """Test that a bytes sequence keeps its reference count."""
        from Bio import cpairwise2

        match_fn = pairwise2.identity_match(1, 0)
        sequence = bytes(_CONVERTIBLE)
        before = sys.getrefcount(sequence)
        cpairwise2._make_score_matrix_fast(
            sequence, _UNCONVERTIBLE, match_fn, -1.0, -0.5, -1.0, -0.5, 0, (1, 1), 1, 1
        )
        self.assertEqual(sys.getrefcount(sequence), before)

    def test_mixed_encodability_does_not_crash(self):
        """Test that aligning one ASCII and one non-ASCII sequence is safe.

        Releasing the converted sequence twice corrupts the allocator rather
        than failing at the point of the mistake, so the interpreter dies
        somewhere unrelated a couple of calls later. Run it in a child process
        so that a regression is reported as a failure here instead of taking
        the whole test run down with it.
        """
        script = (
            "import warnings;warnings.simplefilter('ignore')\n"
            "from Bio import pairwise2\n"
            "for _ in range(10):\n"
            "    pairwise2.align.globalxx('#parotide#gauche#des#',"
            " '#\\u00f9#parotide#gauche#')\n"
        )
        process = subprocess.run(
            [sys.executable, "-c", script],
            capture_output=True,
            text=True,
        )
        self.assertEqual(
            process.returncode,
            0,
            f"child exited with {process.returncode}\n{process.stderr}",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
