# Copyright 2022 by Joao Rodrigues. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Unit tests for the Bio.PDB.CEAligner module."""

import gc
import pickle
import platform
import sys
import sysconfig
import unittest

try:
    import tracemalloc
except ImportError:
    # PyPy has no _tracemalloc, and importing it here would fail the whole
    # module rather than skipping the one test that needs it.
    tracemalloc = None

try:
    import numpy as np
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB."
    ) from None

from Bio.PDB import CEAligner
from Bio.PDB import MMCIFParser
from Bio.PDB import ccealign
from Bio.PDB.ccealign import run_cealign


_stable_cpython_refcounts = (
    platform.python_implementation() == "CPython"
    and not sysconfig.get_config_var("Py_GIL_DISABLED")
)


class CEAlignerTests(unittest.TestCase):
    """Test CEAligner class."""

    @staticmethod
    def _get_ca_coords_as_array(structure):
        xyz_list = [a.coord for a in structure.get_atoms() if a.name == "CA"]
        return np.asarray(xyz_list, dtype=np.float64)

    def test_cealigner(self):
        """Test aligning 7CFN on 6WQA."""
        ref = "PDB/6WQA.cif"
        mob = "PDB/7CFN.cif"
        result = "PDB/7CFN_aligned.cif"

        parser = MMCIFParser(QUIET=1)
        s1 = parser.get_structure("6wqa", ref)
        s2 = parser.get_structure("7cfn", mob)

        aligner = CEAligner()
        aligner.set_reference(s1)
        aligner.align(s2, final_optimization=False)

        self.assertAlmostEqual(aligner.rms, 3.83, places=2)

        # Assert the transformation was done right by comparing
        # the moved coordinates to a 'ground truth' reference.
        # Reference obtained with Pymol's CEAlign code.
        refe = parser.get_structure("7cfn_aligned", result)
        refe_coords = self._get_ca_coords_as_array(refe)
        s2_f_coords = self._get_ca_coords_as_array(s2)

        diff = refe_coords - s2_f_coords
        rmsd = np.sqrt((diff * diff).sum() / len(refe_coords))
        self.assertAlmostEqual(rmsd, 0.0, places=2)

    def test_cealigner_no_transform(self):
        """Test aligning 7CFN on 6WQA without transforming 7CFN."""
        ref = "PDB/6WQA.cif"
        mob = "PDB/7CFN.cif"

        parser = MMCIFParser(QUIET=1)
        s1 = parser.get_structure("6wqa", ref)
        s2 = parser.get_structure("7cfn", mob)

        s2_original_coords = [list(a.coord) for a in s2.get_atoms()]

        aligner = CEAligner()
        aligner.set_reference(s1)
        aligner.align(s2, transform=False, final_optimization=False)
        s2_coords_final = [list(a.coord) for a in s2.get_atoms()]

        self.assertAlmostEqual(aligner.rms, 3.83, places=2)
        self.assertEqual(s2_original_coords, s2_coords_final)

    def test_ce_aligner_final_optimization(self):
        """Test aligning 7CFN on 6WQA with the final optimization."""
        ref = "PDB/6WQA.cif"
        mob = "PDB/7CFN.cif"

        parser = MMCIFParser(QUIET=1)
        s1 = parser.get_structure("6wqa", ref)
        s2 = parser.get_structure("7cfn", mob)

        aligner = CEAligner()
        aligner.set_reference(s1)
        aligner.align(s2)

        self.assertAlmostEqual(aligner.rms, 3.75, places=2)

    def test_cealigner_nucleic(self):
        """Test aligning 1LCD on 1LCD."""
        ref = "PDB/1LCD.cif"
        mob = "PDB/1LCD.cif"

        parser = MMCIFParser(QUIET=1)
        s1 = parser.get_structure("1lcd_ref", ref)
        s2 = parser.get_structure("1lcd_mob", mob)

        aligner = CEAligner()
        aligner.set_reference(s1)
        aligner.align(s2)

        self.assertAlmostEqual(aligner.rms, 0.0, places=3)

    @unittest.skipUnless(
        _stable_cpython_refcounts,
        "GIL-enabled CPython reference counts are required",
    )
    def test_ccealign_reference_ownership(self):
        """Test that run_cealign does not retain stolen references."""
        # What sys.getrefcount reports for a local holding the only reference
        # is not a fixed number: Python 3.14 defers the reference a function
        # local would otherwise own, so it answers 1 where 3.13 answers 2.
        # Measure it from a control object rather than writing it down.
        control = []
        alone = sys.getrefcount(control)

        coords = [[float(i), 0.0, 0.0] for i in range(16)]
        results = run_cealign(coords, coords, 8, 30)
        self.assertEqual(sys.getrefcount(results), alone)

        # Held by the named tuple as well as by the local.
        pair = results[0].path
        self.assertEqual(sys.getrefcount(pair), alone + 1)

        path_a, path_b = pair
        self.assertEqual(sys.getrefcount(path_a), alone + 1)
        self.assertEqual(sys.getrefcount(path_b), alone + 1)

    @unittest.skipUnless(
        _stable_cpython_refcounts and tracemalloc is not None,
        "GIL-enabled CPython memory tracing is required",
    )
    def test_ccealign_path_memory_released(self):
        """Test that run_cealign releases its raw path buffers."""
        coords = [[float(i), 0.0, 0.0] for i in range(40)]
        was_tracing = tracemalloc.is_tracing()
        if not was_tracing:
            tracemalloc.start()
        try:
            for _ in range(5):
                run_cealign(coords, coords, 8, 30)
            gc.collect()
            baseline = tracemalloc.get_traced_memory()[0]

            for _ in range(80):
                run_cealign(coords, coords, 8, 30)
            gc.collect()
            retained = tracemalloc.get_traced_memory()[0] - baseline
        finally:
            if not was_tracing:
                tracemalloc.stop()

        self.assertLess(retained, 128 * 1024)


class RunCEAlignArgumentTests(unittest.TestCase):
    """Argument validation in run_cealign.

    Every case here used to reach the C loops unchecked. The first three
    segfaulted or raised SystemError from inside the extension rather than
    rejecting the input.
    """

    @staticmethod
    def coords(n=40):
        return [[float(i), 0.0, 0.0] for i in range(n)]

    def test_valid_arguments(self):
        """Test that a well formed call still works."""
        coords = self.coords()
        self.assertGreater(len(run_cealign(coords, coords, 8, 30)), 0)

    def test_fragment_size_not_positive(self):
        """Test that a fragment size of zero or less is rejected."""
        coords = self.coords()
        for fragment_size in (0, -1):
            with self.subTest(fragment_size=fragment_size):
                with self.assertRaisesRegex(ValueError, "must be positive"):
                    run_cealign(coords, coords, fragment_size, 30)

    def test_fragment_size_longer_than_input(self):
        """Test that a fragment longer than the structure is rejected."""
        coords = self.coords()
        with self.assertRaisesRegex(ValueError, "at least 100 coordinates"):
            run_cealign(coords, coords, 100, 30)

    def test_coordinate_entry_too_short(self):
        """Test that an entry without three values is rejected."""
        short = [[0.0, 0.0] for _ in range(40)]
        with self.assertRaisesRegex(ValueError, "three values"):
            run_cealign(short, short, 8, 30)

    def test_coordinate_value_not_a_number(self):
        """Test that a non-numeric coordinate is rejected."""
        bad = [["x", 0.0, 0.0] for _ in range(40)]
        with self.assertRaisesRegex(ValueError, "not a number"):
            run_cealign(bad, self.coords(), 8, 30)

    def test_gap_max_out_of_range(self):
        """Test that a gap maximum that cannot be doubled is rejected."""
        coords = self.coords()
        for gap_max in (-1, 2**31 - 1):
            with self.subTest(gap_max=gap_max):
                with self.assertRaisesRegex(ValueError, "maximum gap"):
                    run_cealign(coords, coords, 8, gap_max)

    def test_length_beyond_int_range(self):
        """Test that a reported length that cannot fit an int is rejected.

        The reported length used to be cast straight to a C int, so a
        sequence claiming 2**32 + 40 items reached the C loops as one of
        length 40, bypassing the length checks. On 32-bit builds
        PySequence_Size itself raises OverflowError before the extension's
        own check can, so both exceptions are accepted.
        """

        class LyingLength:
            def __len__(self):
                return 2**32 + 40

            def __getitem__(self, index):
                return (0.0, 0.0, 0.0)

        with self.assertRaises((ValueError, OverflowError)):
            run_cealign(LyingLength(), self.coords(), 8, 30)

    def test_not_a_sequence(self):
        """Test that a non-sequence argument raises TypeError."""
        coords = self.coords()
        with self.assertRaises(TypeError):
            run_cealign(42, coords, 8, 30)

    def test_tuples_are_accepted(self):
        """Test that tuples work as well as lists.

        These used to raise SystemError from PyList_GetItem, because the
        coordinates were read with list-only accessors after a length check
        whose failure was ignored.
        """
        coords = self.coords()
        as_tuples = tuple(tuple(xyz) for xyz in coords)
        self.assertEqual(
            len(run_cealign(as_tuples, as_tuples, 8, 30)),
            len(run_cealign(coords, coords, 8, 30)),
        )


class RunCEAlignResultTypeTests(unittest.TestCase):
    """The CEAlignment result type is shared and picklable.

    run_cealign used to build a fresh CEAlignment type for every result,
    so two results of the same call had distinct types and no result
    could be pickled.
    """

    @staticmethod
    def coords(n=40):
        return [[float(i), 0.0, 0.0] for i in range(n)]

    def test_results_share_one_type(self):
        """Test that every result of a call has the same type."""
        coords = self.coords()
        results = run_cealign(coords, coords, 8, 30)
        self.assertGreater(len(results), 1)
        for result in results[1:]:
            self.assertIs(type(result), type(results[0]))

    def test_result_type_is_module_attribute(self):
        """Test that the result type is ccealign.CEAlignment."""
        coords = self.coords()
        results = run_cealign(coords, coords, 8, 30)
        self.assertIs(type(results[0]), ccealign.CEAlignment)

    def test_result_pickles(self):
        """Test that a result survives a pickle round-trip."""
        coords = self.coords()
        result = run_cealign(coords, coords, 8, 30)[0]
        clone = pickle.loads(pickle.dumps(result))
        self.assertIs(type(clone), type(result))
        self.assertEqual(clone, result)
        self.assertEqual(clone.path, result.path)
        self.assertEqual(clone.z_score, result.z_score)
        self.assertEqual(clone.length, result.length)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
