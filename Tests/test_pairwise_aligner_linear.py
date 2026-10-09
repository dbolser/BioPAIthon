# This code is part of the BioPAIthon distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for the linear-space traceback of PairwiseAligner.align().

Needleman-Wunsch alignments use it where their traceback matrix would be over
512 MiB; the private _set_traceback_limits changes that. It must give the full
traceback matrix's score, number of paths, and paths in the same order.
BIOPAITHON_SLOW_TESTS=1 adds a soak and alignments of over INT_MAX cells.
"""

import contextlib
import os
import platform
import random
import signal
import struct
import sys
import threading
import unittest

try:
    import numpy as np
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align."
    ) from None

from Bio import Align
from Bio import SeqIO
from Bio.Align import _pairwisealigner
from Bio.Align import substitution_matrices
from Bio.Seq import reverse_complement

import support
import test_pairwise_aligner
from memory_growth import requires_growth_measurement

set_limits = _pairwisealigner._set_traceback_limits

# The defaults, as found; _set_traceback_limits returns the previous limits.
DEFAULT = set_limits(None, 1, 1)
set_limits(*DEFAULT)

# Route every Needleman-Wunsch alignment to the linear-space traceback,
# even where it holds more memory than the full matrix.
ALWAYS = -1

SLOW = os.environ.get("BIOPAITHON_SLOW_TESTS") == "1"


@contextlib.contextmanager
def limits(threshold, checkpoint_bytes, block_bytes):
    """Set the traceback limits for the duration of a with block."""
    previous = set_limits(threshold, checkpoint_bytes, block_bytes)
    try:
        yield
    finally:
        set_limits(*previous)


def nw_aligner(**scores):
    """Return a Needleman-Wunsch aligner."""
    aligner = Align.PairwiseAligner(**{"mismatch_score": -1, "gap_score": -1, **scores})
    assert aligner.algorithm == "Needleman-Wunsch", aligner.algorithm
    return aligner


def forced_align(aligner, *args, budgets=(1 << 16, 1 << 12)):
    """Align with the linear-space traceback, and these budgets in bytes."""
    with limits(ALWAYS, *budgets):
        return aligner.align(*args)


def random_dna(rng, length):
    """Return a random DNA sequence as a string."""
    return "".join(rng.choice("ACGT") for _ in range(length))


with limits(ALWAYS, 1, 1):
    LinearPaths = type(nw_aligner().align("A", "A")._paths)
PathGenerator = type(nw_aligner().align("A", "A")._paths)


def coordinates(alignments):
    """Return the coordinates of all alignments, in order, as lists."""
    return [alignment.coordinates.tolist() for alignment in alignments]


def count(alignments):
    """Return len(alignments), or the OverflowError message if it overflows."""
    try:
        return len(alignments)
    except OverflowError as exception:
        return str(exception)


def traced_peak(function, *args):
    """Return the most memory traced while calling function with args."""
    import tracemalloc  # not on PyPy, where the tests using this are skipped

    tracemalloc.start()
    try:
        function(*args)
        return tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()


def random_case(rng, lengths):
    """Return a random Needleman-Wunsch aligner and its align() arguments.

    Draws the alphabet (with or without a wildcard), the scoring (match and
    mismatch scores, an Array with an alphabet, or a NumPy matrix with
    integer-coded sequences), all six gap scores, epsilon and the strand.
    """
    alphabet = rng.choice(["ACGT", "AC", "A"])
    if rng.random() < 0.5:
        alphabet += "N"
    seqA = "".join(rng.choice(alphabet) for _ in range(rng.randint(*lengths)))
    seqB = "".join(rng.choice(alphabet) for _ in range(rng.randint(*lengths)))
    values = [0, 0.5, 1, -0.1, -0.3, -0.7, -1, -2]
    aligner = Align.PairwiseAligner()
    scoring = rng.choice(["compare", "array", "numpy"])
    if scoring == "compare":
        aligner.match_score = rng.choice([1, 2, 1.1, 0.3])
        aligner.mismatch_score = rng.choice([0, -1, -0.3, -0.7])
        if "N" in alphabet:
            aligner.wildcard = "N"
    else:
        size = len(alphabet)
        matrix = np.array([rng.choice(values) for _ in range(size * size)], float)
        matrix = matrix.reshape(size, size)
        if scoring == "array":
            array = substitution_matrices.Array(alphabet, dims=2)
            array[:, :] = matrix
            aligner.substitution_matrix = array
        else:
            aligner.substitution_matrix = matrix
            seqA = np.array([alphabet.index(c) for c in seqA], np.int32)
            seqB = np.array([alphabet.index(c) for c in seqB], np.int32)
    for end in ("left", "internal", "right"):
        for gap in ("insertion", "deletion"):
            setattr(aligner, f"{end}_{gap}_score", rng.choice(values))
    epsilon = rng.choice([None, 0, 0.05, 0.31])
    if epsilon is not None:
        aligner.epsilon = epsilon
    if scoring != "numpy" and rng.random() < 0.5:
        # align() aligns seqA to the reverse complement of the seqB it gets.
        args = (seqA, reverse_complement(seqB), "-")
    else:
        args = (seqA, seqB, "+")
    assert aligner.algorithm == "Needleman-Wunsch", aligner.algorithm
    return aligner, args


def random_budgets(rng, nB):
    """Return a checkpoint and a block budget, from one cell or row to large."""
    row = 8 * (nB + 1)
    checkpoint_bytes = rng.choice([1, row, rng.randint(1, 16 * row), 1 << 30])
    block_bytes = rng.choice([1, nB + 1, rng.randint(1, 16 * (nB + 1)), 1 << 30])
    return checkpoint_bytes, block_bytes


class TestDefaults(unittest.TestCase):
    """Needleman-Wunsch above 512 MiB of traceback matrix, and nothing else."""

    def test_defaults(self):
        # The threshold, the checkpoint budget and the block budget.
        self.assertEqual(DEFAULT, (512 << 20, 32 << 20, 16 << 20))
        self.assertIsNot(LinearPaths, PathGenerator)

    def test_routing(self):
        # The smallest square alignment over the threshold, and one below it.
        pointer = struct.calcsize("P")
        n = 23000
        while (n + 1) * (n + 1 + pointer) <= DEFAULT[0]:
            n += 1
        aligner = nw_aligner()
        rng = random.Random(23166)
        seqA, seqB = random_dna(rng, n), random_dna(rng, n)
        self.assertIs(type(aligner.align(seqA, seqB)._paths), LinearPaths)
        alignments = aligner.align(seqA[:20000], seqB[:2000])
        self.assertIs(type(alignments._paths), PathGenerator)

    def test_only_needleman_wunsch(self):
        # Smith-Waterman, Gotoh, Waterman-Smith-Beyer and FOGSAA keep the
        # full traceback matrix in this step.
        aligners = [Align.PairwiseAligner(mode="local", gap_score=-1)]
        aligners.append(Align.PairwiseAligner(open_gap_score=-2))
        aligners.append(Align.PairwiseAligner(gap_score=lambda i, n: -n))
        aligners.append(
            Align.PairwiseAligner(mode="fogsaa", mismatch_score=-1, gap_score=-2)
        )
        for aligner in aligners:
            with limits(ALWAYS, 1, 1):
                alignments = aligner.align("ACGT", "AGT")
            self.assertIs(type(alignments._paths), PathGenerator, aligner.algorithm)

    def test_threshold(self):
        # The full Needleman-Wunsch matrix of 200 x 200 letters is 201 rows
        # of 201 one-byte cells and a row pointer.
        nbytes = 201 * (201 + struct.calcsize("P"))
        aligner = nw_aligner()
        seqs = ("ACGT" * 50, "AGCT" * 50)
        with limits(nbytes, 1, 1):
            self.assertIs(type(aligner.align(*seqs)._paths), PathGenerator)
        with limits(nbytes - 1, 1, 1):
            self.assertIs(type(aligner.align(*seqs)._paths), LinearPaths)

    def test_set_traceback_limits(self):
        with limits(*DEFAULT):
            self.assertEqual(set_limits(5, 6, 7), DEFAULT)
            self.assertEqual(set_limits(None, 1, 2), (5, 6, 7))
            self.assertEqual(set_limits(ALWAYS, 3, 4), (None, 1, 2))
            self.assertEqual(set_limits(*DEFAULT), (ALWAYS, 3, 4))
            self.assertRaises(ValueError, set_limits, None, 0, 1)
            self.assertRaises(ValueError, set_limits, None, 1, 0)
            self.assertEqual(set_limits(*DEFAULT), DEFAULT)

    def test_laziness(self):
        # repr(), len(), bool() and [0] work in linear space; [1] builds the
        # full matrix.
        aligner = nw_aligner(mismatch_score=0, gap_score=0)
        expected = aligner.align("TACCG", "ACG")
        alignments = forced_align(aligner, "TACCG", "ACG")
        self.assertIs(type(alignments._paths), LinearPaths)
        pointers = (hex(id(expected)), hex(id(alignments)))
        self.assertEqual(repr(alignments), repr(expected).replace(*pointers))
        self.assertEqual(len(alignments), len(expected))
        self.assertTrue(alignments)
        alignments[0]
        self.assertFalse(alignments._paths._materialized)
        alignments[1]
        self.assertTrue(alignments._paths._materialized)
        self.assertEqual(len(alignments), len(expected))


@requires_growth_measurement
class TestMemory(unittest.TestCase):
    """The linear-space traceback is used only where it holds less memory."""

    def peak(self, aligner, seqs, threshold):
        """Return the alignments, and the most memory traced finding [0]."""
        result = []

        def function():
            with limits(threshold, 1, 1):
                alignments = aligner.align(*seqs)
            alignments[0]
            result.append(alignments)

        peak = traced_peak(function)
        return result[0], peak

    def check(self, aligner, seqs, routed):
        """Check routing at threshold 0 against the memory each way takes."""
        alignments, full = self.peak(aligner, seqs, None)
        self.assertIs(type(alignments._paths), PathGenerator)
        alignments, linear = self.peak(aligner, seqs, ALWAYS)
        self.assertIs(type(alignments._paths), LinearPaths)
        shape = (len(seqs[0]), len(seqs[1]))
        self.assertEqual(linear < full, routed, (shape, linear, full))
        with limits(0, 1, 1):
            alignments = aligner.align(*seqs)
        self.assertEqual(type(alignments._paths) is LinearPaths, routed, shape)

    def test_memory(self):
        # A row of trace bits costs less than a row of doubles, so a short
        # target against a long query has too few rows to save.
        aligner = nw_aligner()
        rng = random.Random(8)
        shapes = [(20, 100000, False), (50, 100000, True), (100000, 50, True)]
        for nA, nB, routed in shapes:
            self.check(aligner, (random_dna(rng, nA), random_dna(rng, nB)), routed)

    def test_substitution_matrix(self):
        # The linear-space traceback keeps a copy of the substitution matrix.
        # The full matrix of 1000 x 1000 letters is about 1 MB, and a matrix
        # of 400 x 400 doubles is 1.28 MB.
        rng = np.random.default_rng(400)
        for size, routed in [(100, True), (400, False)]:
            aligner = nw_aligner(substitution_matrix=rng.random((size, size)))
            seqA, seqB = rng.integers(size, size=(2, 1000), dtype=np.int32)
            self.check(aligner, (seqA, seqB), routed)

    def peaks(self, aligner, seqs):
        """Return the most memory traced by align(), repr() and [0], each way."""
        routes = {ALWAYS: LinearPaths, None: PathGenerator}

        def function(threshold):
            with limits(threshold, 1 << 16, 1 << 12):
                alignments = aligner.align(*seqs)
            repr(alignments)
            alignments[0]
            self.assertIs(type(alignments._paths), routes[threshold])

        return [traced_peak(function, threshold) for threshold in routes]

    def test_count(self):
        # Counting the paths for repr() takes a row of counts, not the full
        # matrix of about 9 MB.
        aligner = nw_aligner()
        rng = random.Random(3000)
        seqs = (random_dna(rng, 3000), random_dna(rng, 3000))
        linear, full = self.peaks(aligner, seqs)
        self.assertLess(linear, 2_000_000)
        self.assertGreater(full, 9_000_000)

    def test_skewed(self):
        # Nothing is kept per row of the target but its private copy, of 4
        # bytes a letter, against at least 17 bytes a row for the full
        # matrix.  The arrays are made before the memory is traced.
        aligner = nw_aligner()
        rng = np.random.default_rng(1000000)
        seqA = rng.integers(4, size=1000000, dtype=np.int32)
        seqB = np.array([0, 1, 2, 3, 0, 1, 2, 3], np.int32)
        linear, full = self.peaks(aligner, (seqA, seqB))
        self.assertLess(3 * linear, full)


class TestDifferential(unittest.TestCase):
    """The linear-space traceback against the full traceback matrix."""

    def check(self, aligner, args, budgets):
        message = f"{aligner} {args} {budgets}"
        forced = forced_align(aligner, *args, budgets=budgets)
        expected = aligner.align(*args)
        self.assertIs(type(forced._paths), LinearPaths, message)
        self.assertIs(type(expected._paths), PathGenerator, message)
        self.assertEqual(forced.score, expected.score, message)
        # [0] first, before anything builds the full matrix.
        first = expected[0].coordinates
        self.assertTrue(np.array_equal(forced[0].coordinates, first), message)
        try:
            second = expected[1].coordinates
        except IndexError:
            with self.assertRaises(IndexError):
                forced[1]
        else:
            self.assertTrue(np.array_equal(forced[1].coordinates, second), message)
        # The number of paths, or the same OverflowError message.
        length = count(expected)
        self.assertEqual(count(forced), length, message)
        if isinstance(length, int) and length <= 50:
            self.assertEqual(coordinates(forced), coordinates(expected), message)
        forced.rewind()
        self.assertTrue(np.array_equal(forced[0].coordinates, first), message)
        # len() first, counting in linear space, then [0].
        forced = forced_align(aligner, *args, budgets=budgets)
        self.assertEqual(count(forced), length, message)
        self.assertTrue(np.array_equal(forced[0].coordinates, first), message)
        self.assertFalse(forced._paths._materialized, message)

    def run_cases(self, seed, count, lengths):
        rng = random.Random(seed)
        for _ in range(count):
            aligner, args = random_case(rng, lengths)
            self.check(aligner, args, random_budgets(rng, len(args[1])))

    def test_short(self):
        self.run_cases(seed=2026, count=2000, lengths=(1, 40))

    def test_long(self):
        # Scores that are not dyadic fractions, so that rounding matters.
        rng = random.Random(17)
        values = [-0.1, -0.3, -0.7, -1.3]
        for _ in range(20):
            seqA = random_dna(rng, rng.randint(300, 3000))
            seqB = random_dna(rng, rng.randint(300, 3000))
            aligner = Align.PairwiseAligner(match_score=1.1, mismatch_score=-0.7)
            for end in ("left", "internal", "right"):
                for gap in ("insertion", "deletion"):
                    setattr(aligner, f"{end}_{gap}_score", rng.choice(values))
            aligner.epsilon = rng.choice([1e-6, 0, 0.05, 0.31])
            args = (seqA, seqB, rng.choice("+-"))
            self.check(aligner, args, random_budgets(rng, len(seqB)))

    def test_overflow(self):
        # TestOverflowError of test_pairwise_aligner, in linear space.
        aligner = Align.PairwiseAligner(gap_score=0)
        seqA = SeqIO.read(support.DATA / "Align" / "bsubtilis.fa", "fasta").seq
        seqB = SeqIO.read(support.DATA / "Align" / "ecoli.fa", "fasta").seq
        expected = aligner.align(seqA, seqB)
        forced = forced_align(aligner, seqA, seqB)
        self.assertEqual(
            repr(forced),
            f"<PairwiseAlignments object (>{sys.maxsize} alignments; "
            f"score=1286) at {hex(id(forced))}>",
        )
        self.assertEqual(count(forced), count(expected))
        self.assertRegex(count(forced), "^number of optimal alignments is larger")
        self.assertTrue(np.array_equal(forced[0].coordinates, expected[0].coordinates))
        self.assertFalse(forced._paths._materialized)

    @unittest.skipUnless(SLOW, "set BIOPAITHON_SLOW_TESTS=1 to run")
    def test_soak(self):
        self.run_cases(seed=4589, count=50000, lengths=(1, 60))
        self.run_cases(seed=3407, count=300, lengths=(100, 1500))

    @unittest.skipUnless(SLOW, "set BIOPAITHON_SLOW_TESTS=1 to run")
    def test_more_than_int_max_cells(self):
        # The full matrix would need over 2 GB, so check the path against
        # its own score instead.
        rng = random.Random(47000)
        seqA = random_dna(rng, 47000)
        seqB = list(seqA)
        for _ in range(20):
            position = rng.randrange(len(seqB))
            if rng.random() < 0.5:
                del seqB[position : position + rng.randint(1, 30)]
            else:
                seqB[position:position] = rng.choices("ACGT", k=rng.randint(1, 30))
        seqB = "".join(seqB)
        self.assertGreater((len(seqA) + 1) * (len(seqB) + 1), 2**31 - 1)
        aligner = nw_aligner()
        alignments = forced_align(aligner, seqA, seqB, budgets=DEFAULT[1:])
        alignment = alignments[0]
        counts = alignment.counts()
        score = counts.identities - counts.mismatches - counts.gaps
        self.assertEqual(score, alignments.score)
        self.assertEqual(alignments.score, aligner.score(seqA, seqB))
        self.assertEqual(alignment.coordinates[:, -1].tolist(), [len(seqA), len(seqB)])

    @unittest.skipUnless(SLOW, "set BIOPAITHON_SLOW_TESTS=1 to run")
    def test_read_against_genome(self):
        # A 1,000-letter read against 2.2 million letters, with the default
        # limits: over INT_MAX cells, and a 2.2 GB full matrix.
        rng = random.Random(2200000)
        seqA = random_dna(rng, 2200000)
        seqB = list(seqA[1100000:1101000])
        for _ in range(20):
            seqB[rng.randrange(len(seqB))] = rng.choice("ACGT")
        seqB = "".join(seqB)
        self.assertGreater((len(seqA) + 1) * (len(seqB) + 1), 2**31 - 1)
        aligner = nw_aligner()
        alignments = aligner.align(seqA, seqB)
        self.assertIs(type(alignments._paths), LinearPaths)
        alignment = alignments[0]
        counts = alignment.counts()
        score = counts.identities - counts.mismatches - counts.gaps
        self.assertEqual(score, alignments.score)
        self.assertEqual(alignments.score, aligner.score(seqA, seqB))
        self.assertFalse(alignments._paths._materialized)


class TestSnapshot(unittest.TestCase):
    """Changes made after align() do not change the alignments."""

    def check(self, aligner, seqA, seqB, change):
        expected = aligner.align(seqA, seqB)
        reference = (expected[0].coordinates, expected[1].coordinates, len(expected))
        self.assertGreater(reference[2], 2)
        forced = forced_align(aligner, seqA, seqB, budgets=(1, 1))
        change()
        self.assertTrue(np.array_equal(forced[0].coordinates, reference[0]))
        self.assertTrue(np.array_equal(forced[1].coordinates, reference[1]))
        self.assertEqual(len(forced), reference[2])

    def test_scores(self):
        aligner = nw_aligner(wildcard="N")
        seqA = "GAACTNGACGTTANCC"
        seqB = "GACNTGACGTAACCA"

        def change():
            aligner.match_score = -3
            aligner.mismatch_score = 2
            aligner.gap_score = 0.5
            aligner.epsilon = 10
            aligner.wildcard = "A"

        self.check(aligner, seqA, seqB, change)

    def test_array(self):
        matrix = substitution_matrices.load("NUC.4.4")
        aligner = nw_aligner(substitution_matrix=matrix)
        seqA = "GAACTTGACGTTAGCC"
        seqB = "GACATGACGTAACCA"

        def change():
            for letter in "ACGT":
                matrix[letter, letter] = -10

        self.check(aligner, seqA, seqB, change)

    def test_caller_arrays(self):
        matrix = 2 * np.identity(4) - 1
        aligner = nw_aligner(substitution_matrix=matrix)
        seqA = np.array(["ACGT".index(c) for c in "GAACTTGACGTTAGCCTA"], np.int32)
        seqB = np.array(["ACGT".index(c) for c in "GACATGACGTAACCA"], np.int32)

        def change():
            seqA[:] = 0
            seqB[:] = 1
            matrix[:, :] = 5

        self.check(aligner, seqA, seqB, change)


class TestThreads(unittest.TestCase):
    """Concurrent next(), len() and reset() on one paths object."""

    def test_threads(self):
        aligner = nw_aligner(mismatch_score=0, gap_score=0)
        seqA = "GAACTTGACGTTAGCCTA"
        seqB = "GACATGACGTAACCA"
        expected = aligner.align(seqA, seqB)
        length = len(expected)
        self.assertGreater(length, 10)
        paths = set(expected._paths)
        self.assertEqual(len(paths), length)
        shared = forced_align(aligner, seqA, seqB, budgets=(1, 1))._paths
        barrier = threading.Barrier(4)
        errors = []

        def work(seed):
            rng = random.Random(seed)
            try:
                barrier.wait()
                for _ in range(50):
                    action = rng.choice("nnnlr")
                    if action == "n":
                        try:
                            path = next(shared)
                        except StopIteration:
                            continue
                        if path not in paths:
                            errors.append(f"unexpected path {path}")
                    elif action == "l":
                        if len(shared) != length:
                            errors.append(f"unexpected length {len(shared)}")
                    else:
                        shared.reset()
            except Exception as exception:  # noqa: BLE001
                errors.append(repr(exception))

        threads = [threading.Thread(target=work, args=(seed,)) for seed in range(4)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        self.assertEqual(errors, [])


@unittest.skipUnless(hasattr(signal, "setitimer"), "requires POSIX interval timers")
@unittest.skipUnless(
    platform.python_implementation() == "CPython",
    "PyPy runs Python signal handlers only once a C call has returned",
)
class TestInterrupt(unittest.TestCase):
    """KeyboardInterrupt during the traceback leaves the object usable."""

    @classmethod
    def setUpClass(cls):
        rng = random.Random(3000)
        cls.seqs = (random_dna(rng, 3000), random_dna(rng, 3000))
        cls.aligner = nw_aligner()
        cls.expected = cls.aligner.align(*cls.seqs)

    def interrupt(self, paths, returned, call=next):
        """Call call(paths) with a timer raising KeyboardInterrupt."""
        calls = 0

        def handler(signum, frame):
            nonlocal calls
            calls += 1
            if calls == 3:
                raise KeyboardInterrupt

        previous = signal.signal(signal.SIGALRM, handler)
        signal.setitimer(signal.ITIMER_REAL, 0.001, 0.001)
        try:
            returned.append(call(paths))
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
            signal.signal(signal.SIGALRM, previous)

    def test_interrupt(self):
        expected = self.expected
        reference = (expected[0].coordinates, expected[1].coordinates, count(expected))
        for stage, repeats in (("first path", 20), ("full matrix", 20), ("count", 5)):
            for _ in range(repeats):
                alignments = forced_align(self.aligner, *self.seqs, budgets=DEFAULT[1:])
                paths = alignments._paths
                if stage == "full matrix":
                    next(paths)
                returned = []
                call = len if stage == "count" else next
                with self.assertRaises(KeyboardInterrupt, msg=stage):
                    self.interrupt(paths, returned, call)
                # interrupted inside the C code
                self.assertEqual(returned, [], stage)
                self.assertFalse(paths._materialized, stage)
                alignments.rewind()
                if stage == "count":
                    # counted again in linear space, as nothing was kept
                    self.assertEqual(count(alignments), reference[2])
                    self.assertFalse(paths._materialized)
                self.assertTrue(np.array_equal(alignments[0].coordinates, reference[0]))
                self.assertTrue(np.array_equal(alignments[1].coordinates, reference[1]))
                self.assertEqual(count(alignments), reference[2])

    def test_reentry(self):
        # A signal handler using the alignments while they are being
        # computed gets an error, rather than waiting for itself forever.
        alignments = forced_align(self.aligner, *self.seqs)
        paths = alignments._paths

        def handler(signum, frame):
            len(paths)

        previous = signal.signal(signal.SIGALRM, handler)
        signal.setitimer(signal.ITIMER_REAL, 0.002)
        try:
            with self.assertRaisesRegex(RuntimeError, "being computed"):
                next(paths)
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
            signal.signal(signal.SIGALRM, previous)
        first = self.expected[0].coordinates
        self.assertTrue(np.array_equal(alignments[0].coordinates, first))


# The forced re-run of test_pairwise_aligner.  align() is wrapped: each call
# first checks [0] of a fresh linear-space object against a fresh default one,
# whether or not the test looks at [0], then returns another fresh
# linear-space object.  A call the linear-space traceback does not take is
# made only once, as its gap function may change the sequences.
#
# The Forced classes are made at module level, as pytest ignores load_tests.

original_align = Align.PairwiseAligner.align
align_lock = threading.Lock()  # the tests' threads share the limits
calls = {"tests": 0, "routed": 0}


def checked_align(self, *args, **kwargs):
    """Check [0] against the full matrix, then align in linear space."""
    with align_lock:
        forced = original_align(self, *args, **kwargs)
        if self.algorithm != "Needleman-Wunsch":
            return forced
        # Every Needleman-Wunsch alignment takes the linear-space traceback.
        assert type(forced._paths) is LinearPaths, type(forced._paths)
        calls["routed"] += 1
        with limits(*DEFAULT):
            expected = original_align(self, *args, **kwargs)
        assert forced.score == expected.score, (forced.score, expected.score)
        first, second = forced[0].coordinates, expected[0].coordinates
        assert np.array_equal(first, second), (first, second)
        return original_align(self, *args, **kwargs)


class Forced:
    """Mixin running a test case of test_pairwise_aligner in linear space."""

    def setUp(self):
        calls["tests"] += 1
        previous = set_limits(ALWAYS, 1 << 16, 1 << 12)
        self.addCleanup(set_limits, *previous)
        Align.PairwiseAligner.align = checked_align
        self.addCleanup(setattr, Align.PairwiseAligner, "align", original_align)
        super().setUp()


def forced_classes(module):
    """Return a forced re-run class for each test case class of a module."""
    classes = {}
    for name, value in vars(module).items():
        if (
            isinstance(value, type)
            and issubclass(value, unittest.TestCase)
            and value.__module__ == module.__name__
        ):
            classes[f"Forced{name}"] = type(f"Forced{name}", (Forced, value), {})
    return classes


FORCED = forced_classes(test_pairwise_aligner)
globals().update(FORCED)


def skipped(test):
    """Return whether a decorator skips a test case class or method."""
    return getattr(test, "__unittest_skip__", False)


class TestForcedRerun(unittest.TestCase):
    """Check the forced re-run as a whole; its classes sort before this one."""

    def test_routed(self):
        loader = unittest.TestLoader()
        tests = [
            name
            for cls in FORCED.values()
            if not skipped(cls)
            for name in loader.getTestCaseNames(cls)
            if not skipped(getattr(cls, name))
        ]
        if calls["tests"] < len(tests):
            self.skipTest("needs the whole forced re-run in this process")
        self.assertGreater(calls["routed"], 50)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
