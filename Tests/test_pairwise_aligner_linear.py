# This code is part of the BioPAIthon distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for the linear-space traceback of PairwiseAligner.align().

Needleman-Wunsch, Smith-Waterman and Gotoh global and local alignments use
it where their traceback matrix would be over 512 MiB; the private
_set_traceback_limits changes that.
It must give the full traceback matrix's score, number of paths, and paths in
the same order. BIOPAITHON_SLOW_TESTS=1 adds a soak and alignments of over
INT_MAX cells.
"""

import contextlib
import functools
import itertools
import math
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

import support
import test_pairwise_aligner
from memory_growth import requires_growth_measurement

from Bio import Align
from Bio import SeqIO
from Bio.Align import _pairwisealigner
from Bio.Align import substitution_matrices
from Bio.Seq import reverse_complement

set_limits = _pairwisealigner._set_traceback_limits

# The defaults, as found; _set_traceback_limits returns the previous limits.
DEFAULT = set_limits(None, 1, 1)
set_limits(*DEFAULT)

# Route every Needleman-Wunsch, Smith-Waterman and Gotoh alignment to the
# linear-space traceback, even where it holds more memory than the full
# matrix.
ALWAYS = -1

# No CI job sets it; run the gated tests by hand (about ten minutes).
SLOW = os.environ.get("BIOPAITHON_SLOW_TESTS") == "1"

NW = "Needleman-Wunsch"
GOTOH = "Gotoh global alignment algorithm"
SW = "Smith-Waterman"
GOTOH_LOCAL = "Gotoh local alignment algorithm"
WSB_LOCAL = "Waterman-Smith-Beyer local alignment algorithm"
LOCAL = (SW, GOTOH_LOCAL)
AFFINE = (GOTOH, GOTOH_LOCAL)


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
    assert aligner.algorithm == NW, aligner.algorithm
    return aligner


def gotoh_aligner(**scores):
    """Return a Gotoh global aligner."""
    defaults = {"mismatch_score": -1, "open_gap_score": -2, "extend_gap_score": -1}
    aligner = Align.PairwiseAligner(**{**defaults, **scores})
    assert aligner.algorithm == GOTOH, aligner.algorithm
    return aligner


def sw_aligner(**scores):
    """Return a Smith-Waterman aligner."""
    defaults = {"mode": "local", "mismatch_score": -1, "gap_score": -1}
    aligner = Align.PairwiseAligner(**{**defaults, **scores})
    assert aligner.algorithm == SW, aligner.algorithm
    return aligner


def gotoh_local_aligner(**scores):
    """Return a Gotoh local aligner."""
    defaults = {
        "mode": "local",
        "mismatch_score": -1,
        "open_gap_score": -2,
        "extend_gap_score": -1,
    }
    aligner = Align.PairwiseAligner(**{**defaults, **scores})
    assert aligner.algorithm == GOTOH_LOCAL, aligner.algorithm
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


def coordinates(alignments, limit):
    """Return the coordinates of the first limit alignments, in order, as lists.

    An error ends the list, which then holds its repr in its place.
    """
    result = []
    try:
        for alignment in itertools.islice(alignments, limit):
            result.append(alignment.coordinates.tolist())
    except RuntimeError as exception:
        result.append(repr(exception))
    return result


def outcome(alignments, index):
    """Return the coordinates of alignments[index] as lists, or the error."""
    try:
        return alignments[index].coordinates.tolist()
    except (IndexError, RuntimeError) as exception:
        return repr(exception)


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


def random_case(rng, lengths, algorithm=NW):
    """Return a random aligner and its align() arguments.

    Draws the alphabet (with or without a wildcard), the scoring (match and
    mismatch scores, an Array with an alphabet, or a NumPy matrix with
    integer-coded sequences), all six gap scores (all twelve opening and
    extension scores for Gotoh; the internal ones 0 or less for local
    alignments), epsilon (also half or all of a score, or
    negative, for local alignments) and the strand.
    """
    alphabet = rng.choice(["ACGT", "AC", "A"])
    if rng.random() < 0.5:
        alphabet += "N"
    seqA = "".join(rng.choice(alphabet) for _ in range(rng.randint(*lengths)))
    seqB = "".join(rng.choice(alphabet) for _ in range(rng.randint(*lengths)))
    values = [0, 0.5, 1, -0.1, -0.3, -0.7, -1, -2]
    aligner = Align.PairwiseAligner(mode="local" if algorithm in LOCAL else "global")
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
    while True:
        for end in ("left", "internal", "right"):
            # Local mode rejects a positive internal gap score, and ignores
            # the end gap scores.
            if end == "internal" and algorithm in LOCAL:
                choices = [value for value in values if value <= 0]
            else:
                choices = values
            for gap in ("insertion", "deletion"):
                if algorithm in AFFINE:
                    setattr(aligner, f"open_{end}_{gap}_score", rng.choice(choices))
                    setattr(aligner, f"extend_{end}_{gap}_score", rng.choice(choices))
                else:
                    setattr(aligner, f"{end}_{gap}_score", rng.choice(choices))
        if aligner.algorithm == algorithm:
            break
    epsilons = [None, 0, 0.05, 0.31]
    if algorithm in LOCAL:
        # Scores can then land exactly on a tie's edge.
        epsilons += [0.5, 1, -1e-6, -0.05, -1]
    epsilon = rng.choice(epsilons)
    if epsilon is not None:
        aligner.epsilon = epsilon
    if scoring != "numpy" and rng.random() < 0.5:
        # align() aligns seqA to the reverse complement of the seqB it gets.
        args = (seqA, reverse_complement(seqB), "-")
    else:
        args = (seqA, seqB, "+")
    return aligner, args


def random_budgets(rng, nB, algorithm=NW):
    """Return a checkpoint and a block budget, from one cell or row to large."""
    # Gotoh keeps three scores and two bytes of trace bits a cell.  Local
    # alignments carry six doubles at the start of a saved row (the
    # LinearLocal struct).
    scores, trace_bytes = (3, 2) if algorithm in AFFINE else (1, 1)
    carried = 6 if algorithm in LOCAL else 0
    row = 8 * (carried + scores * (nB + 1))
    trace = trace_bytes * (nB + 1)
    checkpoint_bytes = rng.choice([1, row, rng.randint(1, 16 * row), 1 << 30])
    block_bytes = rng.choice([1, trace, rng.randint(1, 16 * trace), 1 << 30])
    return checkpoint_bytes, block_bytes


def brute_force_sw(seqA, seqB, substitution, insertion, deletion):
    """Return the score and coordinates of the optimal local alignments.

    Enumerates every path, independently of the dynamic programming, for
    tiny sequences: from any cell, by diagonal moves, insertions
    (horizontal) and deletions (vertical), never into row 0 or column 0,
    whose cells are start points only.  A local alignment ends with a
    diagonal move, and every part of it before its end scores above 0 (no
    zero-scoring lead-in).  The score is the highest of any, or 0; the
    optimal alignments score it, and none of them extends a shorter one
    ending with a diagonal move by a part scoring 0 (Waterman & Eggert).
    Use scores that add up exactly in floating point.
    """
    nA = len(seqA)
    nB = len(seqB)
    found = []

    def extend(start, i, j, moves, score, prefix):
        # prefix: the highest score of a part ending with a diagonal move
        for move in "DHV":
            di = move != "H"
            dj = move != "V"
            if i + di > nA or j + dj > nB or i + di == 0 or j + dj == 0:
                continue
            if move == "D":
                step = substitution(seqA[i], seqB[j])
            else:
                step = insertion if move == "H" else deletion
            total = score + step
            if move == "D":
                found.append((total, prefix, start, moves + move))
            if total > 0:
                extend(
                    start,
                    i + di,
                    j + dj,
                    moves + move,
                    total,
                    max(prefix, total) if move == "D" else prefix,
                )

    for i in range(nA + 1):
        for j in range(nB + 1):
            extend((i, j), i, j, "", 0, -math.inf)
    best = max([total for total, prefix, start, moves in found if total > 0] or [0])
    paths = {
        path_coordinates(start, moves)
        for total, prefix, start, moves in found
        if best > 0 and total == best and prefix < best
    }
    return best, paths


def path_coordinates(start, moves):
    """Return the coordinates of a path from start by moves, as tuples.

    Moves are "D" (diagonal), "H" (an insertion) and "V" (a deletion).
    """
    i, j = start
    points = [(i, j)]
    for k, move in enumerate(moves):
        i += move != "H"
        j += move != "V"
        if k > 0 and move == moves[k - 1]:
            points[-1] = (i, j)
        else:
            points.append((i, j))
    return tuple(zip(*points))


def brute_force_gotoh_local(
    seqA, seqB, substitution, open_A, extend_A, open_B, extend_B
):
    """Return the score and coordinates of the optimal Gotoh local alignments.

    As brute_force_sw, with affine gap scores: a gap opens with the opening
    score, also after a gap in the other sequence, and goes on with the
    extension score (A for insertions, B for deletions).  As in
    GOTOH_LOCAL_ALIGN, a path starts in any cell, and makes no gap into
    row 0, column 0, the last row or the last column.  The parts of a path
    before its end score above 0, and the optimal alignments extend no
    shorter one ending with a diagonal move by a part scoring 0.
    """
    nA = len(seqA)
    nB = len(seqB)
    found = []

    def extend(start, i, j, moves, score, prefix):
        for move in "DHV":
            ni = i + (move != "H")
            nj = j + (move != "V")
            if ni > nA or nj > nB:
                continue
            if move == "D":
                step = substitution(seqA[i], seqB[j])
            elif ni == 0 or nj == 0 or ni == nA or nj == nB:
                continue
            elif move == "H":
                step = extend_A if moves.endswith("H") else open_A
            else:
                step = extend_B if moves.endswith("V") else open_B
            total = score + step
            if move == "D":
                found.append((total, prefix, start, moves + move))
            if total > 0:
                extend(
                    start,
                    ni,
                    nj,
                    moves + move,
                    total,
                    max(prefix, total) if move == "D" else prefix,
                )

    for i in range(nA + 1):
        for j in range(nB + 1):
            extend((i, j), i, j, "", 0, -math.inf)
    best = max([total for total, prefix, start, moves in found if total > 0] or [0])
    paths = {
        path_coordinates(start, moves)
        for total, prefix, start, moves in found
        if best > 0 and total == best and prefix < best
    }
    return best, paths


def tied_gotoh_local(
    seqA, seqB, substitution, open_A, extend_A, open_B, extend_B, epsilon
):
    """Return the score and coordinates of Gotoh local alignments with ties.

    With an epsilon as large as the steps between scores, a score in M, Ix
    or Iy is that of the first of the matrices before it, in that order,
    that no other beats by more than epsilon, and its trace is all of them
    less than epsilon below that; M below epsilon is a start point, and Ix
    and Iy below epsilon are no path.  The score is the highest of M, or 0;
    the end points are the cells of M less than epsilon below it, and no
    path goes on from one.  For a positive epsilon below the steps between
    scores, it agrees with brute_force_gotoh_local.
    """
    nA = len(seqA)
    nB = len(seqB)
    scores = {}
    traces = {}

    def select(*candidates):
        (score, matrix), *others = candidates
        trace = {matrix}
        for value, matrix in others:
            if value > score + epsilon:
                score, trace = value, {matrix}
            elif value > score - epsilon:
                trace.add(matrix)
        return score, trace

    for i in range(nA + 1):
        for j in range(nB + 1):
            if i == 0 or j == 0:
                scores[i, j] = {"M": 0, "Ix": -math.inf, "Iy": -math.inf}
                traces[i, j] = {"M": None}
                continue
            last = scores[i - 1, j - 1]
            score, trace = select(*((last[m], m) for m in ("M", "Ix", "Iy")))
            score += substitution(seqA[i - 1], seqB[j - 1])
            if score < epsilon:
                score, trace = 0, None  # a start point
            scores[i, j] = {"M": score, "Ix": -math.inf, "Iy": -math.inf}
            traces[i, j] = {"M": trace, "Ix": set(), "Iy": set()}
            if i == nA or j == nB:
                continue
            for gap, last, opening, extension in (
                ("Ix", scores[i - 1, j], open_B, extend_B),
                ("Iy", scores[i, j - 1], open_A, extend_A),
            ):
                score, trace = select(
                    *(
                        (last[m] + (extension if m == gap else opening), m)
                        for m in ("M", "Ix", "Iy")
                    )
                )
                if score >= epsilon:
                    scores[i, j][gap] = score
                    traces[i, j][gap] = trace
    starts = [cell for cell in traces if traces[cell]["M"] is None]
    best = max([s["M"] for cell, s in scores.items() if cell not in starts] or [0])
    ends = {
        cell
        for cell, s in scores.items()
        if cell not in starts and s["M"] > best - epsilon
    }
    moves = {"M": "D", "Ix": "V", "Iy": "H"}
    paths = set()

    def trace_back(i, j, matrix, path):
        trace = traces[i, j][matrix]
        if trace is None:
            paths.add(path_coordinates((i, j), path))
            return
        move = moves[matrix]
        i -= move != "H"
        j -= move != "V"
        for matrix in trace:
            if matrix != "M" or (i, j) not in ends:
                trace_back(i, j, matrix, move + path)

    for i, j in ends:
        trace_back(i, j, "M", "")
    return best, paths


def brute_force_wsb_local(seqA, seqB, substitution, insertion, deletion):
    """Return the score and coordinates of the optimal WSB local alignments.

    As brute_force_gotoh_local, with gap score functions of the length of
    the gap (insertion for gaps in A, deletion for gaps in B): a gap is one
    block, scoring the function of its length, and follows a diagonal move
    or a gap in the other sequence.  The parts of a path before its end,
    ending with a diagonal move or a whole gap, score above 0.
    """
    nA = len(seqA)
    nB = len(seqB)
    found = []

    def extend(start, i, j, moves, score, prefix):
        if i < nA and j < nB:
            total = score + substitution(seqA[i], seqB[j])
            found.append((total, prefix, start, moves + "D"))
            if total > 0:
                extend(start, i + 1, j + 1, moves + "D", total, max(prefix, total))
        for move, function in (("V", deletion), ("H", insertion)):
            if moves.endswith(move):
                continue
            for k in range(1, max(nA, nB)):
                ni = i + k * (move == "V")
                nj = j + k * (move == "H")
                # no gap into row 0, column 0, the last row or the last column
                if not (0 < ni < nA and 0 < nj < nB):
                    break
                total = score + function(k)
                if total > 0:
                    extend(start, ni, nj, moves + move * k, total, prefix)

    for i in range(nA + 1):
        for j in range(nB + 1):
            extend((i, j), i, j, "", 0, -math.inf)
    best = max([total for total, prefix, start, moves in found if total > 0] or [0])
    paths = {
        path_coordinates(start, moves)
        for total, prefix, start, moves in found
        if best > 0 and total == best and prefix < best
    }
    return best, paths


def tied_wsb_local(seqA, seqB, substitution, insertion, deletion, epsilon):
    """Return the score and coordinates of WSB local alignments with ties.

    As tied_gotoh_local, with gap score functions of the length of the gap:
    a score in Ix (or Iy) is that of the first of M and Iy (or Ix) of each
    gap length in turn that no other beats by more than epsilon.  For a
    positive epsilon below the steps between scores, it agrees with
    brute_force_wsb_local.
    """
    nA = len(seqA)
    nB = len(seqB)
    scores = {}
    traces = {}

    def select(*candidates):
        (score, source), *others = candidates
        trace = {source}
        for value, source in others:
            if value > score + epsilon:
                score, trace = value, {source}
            elif value > score - epsilon:
                trace.add(source)
        return score, trace

    for i in range(nA + 1):
        for j in range(nB + 1):
            if i == 0 or j == 0:
                scores[i, j] = {"M": 0, "Ix": -math.inf, "Iy": -math.inf}
                traces[i, j] = {"M": None}
                continue
            last = scores[i - 1, j - 1]
            score, trace = select(*((last[m], (m, 1)) for m in ("M", "Ix", "Iy")))
            score += substitution(seqA[i - 1], seqB[j - 1])
            if score < epsilon:
                score, trace = 0, None  # a start point
            scores[i, j] = {"M": score, "Ix": -math.inf, "Iy": -math.inf}
            traces[i, j] = {"M": trace, "Ix": set(), "Iy": set()}
            if i == nA or j == nB:
                continue
            for gap, other, length, function in (
                ("Ix", "Iy", i, deletion),
                ("Iy", "Ix", j, insertion),
            ):
                candidates = []
                for k in range(1, length + 1):
                    cell = (i - k, j) if gap == "Ix" else (i, j - k)
                    for m in ("M", other):
                        candidates.append((scores[cell][m] + function(k), (m, k)))
                score, trace = select(*candidates)
                if score >= epsilon:
                    scores[i, j][gap] = score
                    traces[i, j][gap] = trace
    starts = [cell for cell in traces if traces[cell]["M"] is None]
    best = max([s["M"] for cell, s in scores.items() if cell not in starts] or [0])
    ends = {
        cell
        for cell, s in scores.items()
        if cell not in starts and s["M"] > best - epsilon
    }
    moves = {"M": "D", "Ix": "V", "Iy": "H"}
    paths = set()

    def trace_back(i, j, matrix, path):
        trace = traces[i, j][matrix]
        if trace is None:
            paths.add(path_coordinates((i, j), path))
            return
        move = moves[matrix]
        for source, k in trace:
            ni = i - k * (move != "H")
            nj = j - k * (move != "V")
            if source != "M" or (ni, nj) not in ends:
                trace_back(ni, nj, source, move * k + path)

    for i, j in ends:
        trace_back(i, j, "M", "")
    return best, paths


class TestDefaults(unittest.TestCase):
    """Alignments above 512 MiB of traceback matrix, and nothing else."""

    def test_defaults(self):
        # The threshold, the checkpoint budget and the block budget.
        self.assertEqual(DEFAULT, (512 << 20, 32 << 20, 16 << 20))
        self.assertIsNot(LinearPaths, PathGenerator)

    def test_routing(self):
        # The smallest square alignment over the threshold, and one below it,
        # with the default budgets: outside BIOPAITHON_SLOW_TESTS, the only
        # check of the linear-space traceback's memory estimate against a
        # 512 MiB matrix.  A Needleman-Wunsch or Smith-Waterman matrix row
        # has a byte a cell and a row pointer, a Gotoh one two of each.  One
        # kernel of each size: test_threshold pins the sizing of all four,
        # and a local kernel's estimate differs only by six doubles a row.
        pointer = struct.calcsize("P")
        for aligner, size, n in (
            (nw_aligner(), 1, 23000),
            (gotoh_aligner(), 2, 16000),
        ):
            while (n + 1) * size * (n + 1 + pointer) <= DEFAULT[0]:
                n += 1
            rng = random.Random(n)
            seqA, seqB = random_dna(rng, n), random_dna(rng, n)
            alignments = aligner.align(seqA, seqB)
            self.assertIs(type(alignments._paths), LinearPaths, aligner.algorithm)
            alignments = aligner.align(seqA[:20000], seqB[:2000])
            self.assertIs(type(alignments._paths), PathGenerator, aligner.algorithm)

    def test_algorithms(self):
        # Needleman-Wunsch, Smith-Waterman and Gotoh alignments can be
        # routed; Waterman-Smith-Beyer and FOGSAA keep the full traceback
        # matrix.
        aligners = [
            (make(), LinearPaths)
            for make in (nw_aligner, gotoh_aligner, sw_aligner, gotoh_local_aligner)
        ]
        for aligner in (
            Align.PairwiseAligner(mode="local", gap_score=lambda i, n: -n),
            Align.PairwiseAligner(gap_score=lambda i, n: -n),
            Align.PairwiseAligner(mode="fogsaa", mismatch_score=-1, gap_score=-2),
        ):
            aligners.append((aligner, PathGenerator))
        for aligner, expected in aligners:
            with limits(ALWAYS, 1, 1):
                alignments = aligner.align("ACGT", "AGT")
            self.assertIs(type(alignments._paths), expected, aligner.algorithm)

    def test_threshold(self):
        # The full Needleman-Wunsch or Smith-Waterman matrix of 200 x 200
        # letters is 201 rows of 201 one-byte cells and a row pointer; the
        # Gotoh one has 201 two-byte cells and two row pointers a row.
        pointer = struct.calcsize("P")
        seqs = ("ACGT" * 50, "AGCT" * 50)
        for aligner, size in (
            (nw_aligner(), 1),
            (gotoh_aligner(), 2),
            (sw_aligner(), 1),
            (gotoh_local_aligner(), 2),
        ):
            nbytes = 201 * size * (201 + pointer)
            with limits(nbytes, 1, 1):
                alignments = aligner.align(*seqs)
            self.assertIs(type(alignments._paths), PathGenerator, aligner.algorithm)
            with limits(nbytes - 1, 1, 1):
                alignments = aligner.align(*seqs)
            self.assertIs(type(alignments._paths), LinearPaths, aligner.algorithm)

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
        for aligner in (
            nw_aligner(mismatch_score=0, gap_score=0),
            gotoh_aligner(mismatch_score=0, open_gap_score=-0.5, extend_gap_score=0),
            sw_aligner(mismatch_score=0, gap_score=0),
            gotoh_local_aligner(
                mismatch_score=0, open_gap_score=-0.5, extend_gap_score=0
            ),
        ):
            expected = aligner.align("TACCG", "ACG")
            self.assertGreater(len(expected), 1)
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
        # target against a long query has too few rows to save.  Gotoh keeps
        # three doubles a column for two bytes of trace bits.  (Between
        # these shapes are close calls, which stay on the faster full
        # matrix.)
        rng = random.Random(8)
        shapes = [(20, 100000, False), (50, 100000, True), (100000, 50, True)]
        for nA, nB, routed in shapes:
            seqs = (random_dna(rng, nA), random_dna(rng, nB))
            self.check(nw_aligner(), seqs, routed)
        shapes = [(20, 100000, False), (100, 100000, True), (100000, 50, True)]
        for nA, nB, routed in shapes:
            seqs = (random_dna(rng, nA), random_dna(rng, nB))
            self.check(gotoh_aligner(), seqs, routed)
        shapes = [(20, 100000, False), (50, 100000, True), (100000, 50, True)]
        for nA, nB, routed in shapes:
            seqs = (random_dna(rng, nA), random_dna(rng, nB))
            self.check(sw_aligner(), seqs, routed)
        shapes = [(20, 100000, False), (100, 100000, True), (100000, 50, True)]
        for nA, nB, routed in shapes:
            seqs = (random_dna(rng, nA), random_dna(rng, nB))
            self.check(gotoh_local_aligner(), seqs, routed)

    def test_substitution_matrix(self):
        # The linear-space traceback keeps a copy of the substitution matrix.
        # The full matrix of 1000 x 1000 letters is about 1 MB, or 2 MB for
        # Gotoh, and a matrix of 400 x 400 doubles is 1.28 MB, of 600 x 600
        # doubles 2.88 MB.
        rng = np.random.default_rng(400)
        for make, size, routed in [
            (nw_aligner, 100, True),
            (nw_aligner, 400, False),
            (gotoh_aligner, 100, True),
            (gotoh_aligner, 600, False),
            (sw_aligner, 100, True),
            (sw_aligner, 400, False),
            (gotoh_local_aligner, 100, True),
            (gotoh_local_aligner, 600, False),
        ]:
            aligner = make(substitution_matrix=rng.random((size, size)))
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
        # Counting the paths for repr() takes a row of counts (three for
        # Gotoh), not the full matrix of about 9 MB (18 MB for Gotoh).
        rng = random.Random(3000)
        seqs = (random_dna(rng, 3000), random_dna(rng, 3000))
        for aligner, size in (
            (nw_aligner(), 9_000_000),
            (gotoh_aligner(), 18_000_000),
            (sw_aligner(), 9_000_000),
            (gotoh_local_aligner(), 18_000_000),
        ):
            linear, full = self.peaks(aligner, seqs)
            self.assertLess(linear, 2_000_000, aligner.algorithm)
            self.assertGreater(full, size, aligner.algorithm)

    def test_skewed(self):
        # Nothing is kept per row of the target but its private copy, of 4
        # bytes a letter, against at least 17 bytes a row for the full
        # matrix (34 for Gotoh).  The arrays are made before the memory is
        # traced.
        rng = np.random.default_rng(1000000)
        seqA = rng.integers(4, size=1000000, dtype=np.int32)
        seqB = np.array([0, 1, 2, 3, 0, 1, 2, 3], np.int32)
        for aligner in (
            nw_aligner(),
            gotoh_aligner(),
            sw_aligner(),
            gotoh_local_aligner(),
        ):
            linear, full = self.peaks(aligner, (seqA, seqB))
            self.assertLess(3 * linear, full, aligner.algorithm)


class TestDifferential(unittest.TestCase):
    """The linear-space traceback against the full traceback matrix."""

    def check(self, aligner, args, budgets):
        message = f"{aligner} {args} {budgets}"
        forced = forced_align(aligner, *args, budgets=budgets)
        expected = aligner.align(*args)
        self.assertIs(type(forced._paths), LinearPaths, message)
        self.assertIs(type(expected._paths), PathGenerator, message)
        self.assertEqual(forced.score, expected.score, message)
        # [0] first, before anything builds the full matrix; then [1].  Each
        # is the same alignment, or the same error: no alignment, or one of
        # the full matrix's PathGenerator.
        first = outcome(expected, 0)
        self.assertEqual(outcome(forced, 0), first, message)
        self.assertEqual(outcome(forced, 1), outcome(expected, 1), message)
        # The number of paths, or the same OverflowError message.
        length = count(expected)
        self.assertEqual(count(forced), length, message)
        if isinstance(length, int) and length <= 50:
            # One more than len(), to see that neither goes on past it.
            limit = length + 1
            paths = coordinates(expected, limit)
            self.assertEqual(coordinates(forced, limit), paths, message)
            self.assertEqual(len(paths), length, message)
        # The same after a reset (see test_gotoh_local_score_epsilon).
        forced.rewind()
        expected.rewind()
        self.assertEqual(outcome(expected, 0), first, message)
        self.assertEqual(outcome(forced, 0), first, message)
        # len() first, counting in linear space, then [0].
        forced = forced_align(aligner, *args, budgets=budgets)
        self.assertEqual(count(forced), length, message)
        self.assertEqual(outcome(forced, 0), first, message)
        self.assertFalse(forced._paths._materialized, message)

    def run_cases(self, seed, count, lengths, algorithm=NW):
        rng = random.Random(seed)
        for _ in range(count):
            aligner, args = random_case(rng, lengths, algorithm)
            self.check(aligner, args, random_budgets(rng, len(args[1]), algorithm))

    def test_short(self):
        self.run_cases(seed=2026, count=2000, lengths=(1, 40))

    def test_short_gotoh(self):
        self.run_cases(seed=2027, count=2000, lengths=(1, 40), algorithm=GOTOH)

    def test_short_sw(self):
        self.run_cases(seed=2028, count=2000, lengths=(1, 40), algorithm=SW)

    def test_short_gotoh_local(self):
        self.run_cases(seed=2029, count=2000, lengths=(1, 40), algorithm=GOTOH_LOCAL)

    def test_sw_brute_force(self):
        # The full matrix's score, number of paths and paths against every
        # optimal local alignment of tiny sequences, found by brute force,
        # with negative and zero gap scores (local mode rejects positive
        # ones); then the linear-space traceback against the full matrix.
        rng = random.Random(2030)
        values = [-2, -1, -0.5, 0, 0.5, 1]
        gap_values = [value for value in values if value <= 0]
        for _ in range(1500):
            alphabet = rng.choice(["AC", "ACG"])
            seqA = "".join(rng.choice(alphabet) for _ in range(rng.randint(1, 5)))
            seqB = "".join(rng.choice(alphabet) for _ in range(rng.randint(1, 5)))
            insertion = rng.choice(gap_values)
            deletion = rng.choice(gap_values) if rng.random() < 0.3 else insertion
            if rng.random() < 0.5:
                match = rng.choice([1, 2, 0.5])
                mismatch = rng.choice([-1, 0, -0.5])
                aligner = sw_aligner(match_score=match, mismatch_score=mismatch)

                def substitution(a, b, match=match, mismatch=mismatch):
                    return match if a == b else mismatch

            else:
                matrix = substitution_matrices.Array(alphabet, dims=2)
                matrix[:, :] = [
                    [rng.choice(values) for _ in alphabet] for _ in alphabet
                ]
                aligner = sw_aligner(substitution_matrix=matrix)

                def substitution(a, b, matrix=matrix):
                    return matrix[a, b]

            aligner.insertion_score = insertion
            aligner.deletion_score = deletion
            # As exact, as the scores are dyadic and the edges of ties
            # further apart than epsilon.
            aligner.epsilon = rng.choice([1e-6, 0.25])
            message = f"{aligner} {seqA} {seqB}"
            score, paths = brute_force_sw(seqA, seqB, substitution, insertion, deletion)
            alignments = aligner.align(seqA, seqB)
            self.assertEqual(alignments.score, score, message)
            self.assertEqual(aligner.score(seqA, seqB), score, message)
            self.assertEqual(len(alignments), len(paths), message)
            found = [
                tuple(map(tuple, coordinates))
                for coordinates in coordinates(alignments, len(paths) + 1)
            ]
            self.assertEqual(len(found), len(paths), message)
            self.assertEqual(set(found), paths, message)
            self.check(aligner, (seqA, seqB, "+"), random_budgets(rng, len(seqB), SW))

    def test_gotoh_local_brute_force(self):
        # As test_sw_brute_force, for Gotoh local alignments, with epsilon
        # from below the steps between scores, where brute force finds the
        # paths, to far above, where they follow the ties (tied_gotoh_local).
        # score() ignores epsilon, so it always agrees with brute force.
        rng = random.Random(2031)
        values = [-2, -1, -0.5, 0, 0.5, 1]
        # local mode rejects positive gap scores
        gap_values = [value for value in values if value <= 0]
        for _ in range(1500):
            alphabet = rng.choice(["AC", "ACG"])
            seqA = "".join(rng.choice(alphabet) for _ in range(rng.randint(1, 5)))
            seqB = "".join(rng.choice(alphabet) for _ in range(rng.randint(1, 5)))
            if rng.random() < 0.5:
                match = rng.choice([1, 2, 0.5])
                mismatch = rng.choice([-1, 0, -0.5])
                scores = {"match_score": match, "mismatch_score": mismatch}

                def substitution(a, b, match=match, mismatch=mismatch):
                    return match if a == b else mismatch

            else:
                matrix = substitution_matrices.Array(alphabet, dims=2)
                matrix[:, :] = [
                    [rng.choice(values) for _ in alphabet] for _ in alphabet
                ]
                scores = {"substitution_matrix": matrix}

                def substitution(a, b, matrix=matrix):
                    return matrix[a, b]

            while True:
                gaps = [rng.choice(gap_values) for _ in range(4)]
                if rng.random() < 0.7:
                    gaps[2:] = gaps[:2]
                open_A, extend_A, open_B, extend_B = gaps
                aligner = Align.PairwiseAligner(
                    mode="local",
                    open_insertion_score=open_A,
                    extend_insertion_score=extend_A,
                    open_deletion_score=open_B,
                    extend_deletion_score=extend_B,
                    **scores,
                )
                if aligner.algorithm == GOTOH_LOCAL:
                    break
            # The scores are dyadic, so they add up exactly.
            aligner.epsilon = rng.choice([1e-6, 0.25, 0.5, 1, 2, 5])
            message = f"{aligner} {seqA} {seqB}"
            best, paths = brute_force_gotoh_local(seqA, seqB, substitution, *gaps)
            self.assertEqual(aligner.score(seqA, seqB), best, message)
            tied = tied_gotoh_local(seqA, seqB, substitution, *gaps, aligner.epsilon)
            if aligner.epsilon < 0.5:
                self.assertEqual(tied, (best, paths), message)
            best, paths = tied
            alignments = aligner.align(seqA, seqB)
            self.assertEqual(alignments.score, best, message)
            self.assertEqual(len(alignments), len(paths), message)
            found = [
                tuple(map(tuple, coordinates))
                for coordinates in coordinates(alignments, len(paths) + 1)
            ]
            self.assertEqual(len(found), len(paths), message)
            self.assertEqual(set(found), paths, message)
            self.check(
                aligner,
                (seqA, seqB, "+"),
                random_budgets(rng, len(seqB), GOTOH_LOCAL),
            )

    def run_long_cases(self, seed, count, algorithm):
        # Scores that are not dyadic fractions, so that rounding matters.
        rng = random.Random(seed)
        values = [-0.1, -0.3, -0.7, -1.3]
        mode = "local" if algorithm in LOCAL else "global"
        for _ in range(count):
            seqA = random_dna(rng, rng.randint(300, 3000))
            seqB = random_dna(rng, rng.randint(300, 3000))
            aligner = Align.PairwiseAligner(
                mode=mode, match_score=1.1, mismatch_score=-0.7
            )
            for end in ("left", "internal", "right"):
                for gap in ("insertion", "deletion"):
                    if algorithm in AFFINE:
                        setattr(aligner, f"open_{end}_{gap}_score", rng.choice(values))
                        setattr(
                            aligner, f"extend_{end}_{gap}_score", rng.choice(values)
                        )
                    else:
                        setattr(aligner, f"{end}_{gap}_score", rng.choice(values))
            if aligner.algorithm != algorithm:
                aligner.open_internal_insertion_score = -2.9
            aligner.epsilon = rng.choice([1e-6, 0, 0.05, 0.31])
            args = (seqA, seqB, rng.choice("+-"))
            self.check(aligner, args, random_budgets(rng, len(seqB), algorithm))

    def test_long(self):
        self.run_long_cases(seed=17, count=20, algorithm=NW)

    def test_long_gotoh(self):
        self.run_long_cases(seed=18, count=12, algorithm=GOTOH)

    def test_long_sw(self):
        self.run_long_cases(seed=19, count=12, algorithm=SW)

    def test_long_gotoh_local(self):
        self.run_long_cases(seed=20, count=12, algorithm=GOTOH_LOCAL)

    def test_read(self):
        # A read with a few changes, inside a longer random target, as a
        # local alignment is mostly used: the path ends inside the matrix,
        # and starts well away from row 0 and column 0.
        rng = random.Random(1912)
        for _ in range(8):
            seqA = random_dna(rng, rng.randint(1000, 4000))
            start = rng.randrange(len(seqA) - 500)
            seqB = list(seqA[start : start + rng.randint(100, 500)])
            for _ in range(len(seqB) // 20):
                position = rng.randrange(len(seqB))
                change = rng.choice("sid")
                if change == "s":
                    seqB[position] = rng.choice("ACGT")
                elif change == "i":
                    seqB.insert(position, rng.choice("ACGT"))
                else:
                    del seqB[position]
            seqB = random_dna(rng, 50) + "".join(seqB) + random_dna(rng, 50)
            algorithm = rng.choice(LOCAL)
            make = sw_aligner if algorithm == SW else gotoh_local_aligner
            aligner = make(match_score=rng.choice([1, 2, 1.1]))
            aligner.epsilon = rng.choice([1e-6, 0, 0.05])
            args = (seqA, seqB, "+")
            if rng.random() < 0.5:
                args = (seqA, reverse_complement(seqB), "-")
            self.check(aligner, args, random_budgets(rng, len(seqB), algorithm))

    def test_presets(self):
        # The scoring schemes using Gotoh: blastn on DNA, blastp on protein,
        # globally and locally.
        rng = random.Random(1255)
        for scoring, alphabet, strands in (
            ("blastn", "ACGT", "+-"),
            ("blastp", "ACDEFGHIKLMNPQRSTVWY", "+"),
        ):
            for mode, algorithm in (("global", GOTOH), ("local", GOTOH_LOCAL)):
                aligner = Align.PairwiseAligner(scoring=scoring, mode=mode)
                self.assertEqual(aligner.algorithm, algorithm)
                for _ in range(4):
                    seqA = "".join(rng.choices(alphabet, k=rng.randint(100, 1500)))
                    seqB = "".join(rng.choices(alphabet, k=rng.randint(100, 1500)))
                    args = (seqA, seqB, rng.choice(strands))
                    budgets = random_budgets(rng, len(seqB), algorithm)
                    self.check(aligner, args, budgets)

    def test_stopped(self):
        # Infinite or NaN scores, or an infinite epsilon, leave traces into
        # column 0 or row 0, or the origin's Ix or Iy, where no path
        # starts.  GOTOH_GLOBAL_ALIGN's clean-up removes them, but the walk
        # follows them and stops, and asks the full matrix, whose path goes
        # elsewhere, if anywhere.  (The full matrix's path used to stop
        # there too, a tuple holding NULLs that crashed Python, while len()
        # did not count it.)
        no_end_deletions = gotoh_aligner(end_deletion_score=-math.inf)
        no_end_insertions = gotoh_aligner(end_insertion_score=-math.inf)
        no_gaps = gotoh_aligner(open_gap_score=-math.inf)
        nan = gotoh_aligner(open_gap_score=math.nan)
        infinite_gap = gotoh_aligner(
            gap_score=math.inf, open_left_gap_score=-1, epsilon=-math.inf
        )
        every_tie = gotoh_aligner(epsilon=math.inf)
        for aligner, seqA, seqB, paths in [
            # stops in column 0, in a block
            (no_end_deletions, "AA", "A", 0),
            (no_gaps, "AACGTTGCA", "ACG", 0),
            (nan, "GGGGAACGT", "AACG", 0),
            # stops in row 0
            (no_end_insertions, "A", "AA", 0),
            (no_gaps, "ACG", "AACGTTGCA", 0),
            (nan, "AACG", "GGGGAACGT", 0),
            # stops at the origin, in Iy
            (infinite_gap, "GT", "A", 0),
            # stops in column 0, where the full matrix's path goes elsewhere
            (every_tie, "ACGT", "AGT", 129),
        ]:
            message = f"{aligner.algorithm} {seqA} {seqB}"
            expected = aligner.align(seqA, seqB)
            forced = forced_align(aligner, seqA, seqB, budgets=(1, 1))
            self.assertIs(type(forced._paths), LinearPaths, message)
            self.assertEqual(repr(forced.score), repr(expected.score), message)
            self.assertEqual(len(expected), paths, message)
            self.assertEqual(len(forced), paths, message)
            self.assertFalse(forced._paths._materialized, message)
            found = [alignment.coordinates.tolist() for alignment in expected]
            self.assertEqual(len(found), paths, message)
            self.assertEqual(outcome(forced, 0), outcome(expected, 0), message)
            self.assertTrue(forced._paths._materialized, message)
            self.assertEqual(
                [alignment.coordinates.tolist() for alignment in forced],
                found,
                message,
            )

    def test_no_paths(self):
        # A negative epsilon clears every trace of Gotoh's last cell, so the
        # full matrix has no paths, and nor has the linear-space traceback.
        # Needleman-Wunsch keeps its last trace, and so its paths.
        for aligner, seqA, seqB in [
            (gotoh_aligner(epsilon=-1), "ACGT", "AGT"),
            (gotoh_aligner(epsilon=-1e-6), "GAACT", "GAT"),
            (Align.PairwiseAligner(scoring="blastn", epsilon=-1e-9), "GAACT", "GAT"),
        ]:
            message = f"{aligner.algorithm} {seqA} {seqB} {aligner.epsilon}"
            expected = aligner.align(seqA, seqB)
            forced = forced_align(aligner, seqA, seqB, budgets=(1, 1))
            self.assertIs(type(forced._paths), LinearPaths, message)
            self.assertEqual(forced.score, expected.score, message)
            self.assertEqual(len(expected), 0, message)
            self.assertEqual(list(expected), [], message)
            with self.assertRaises(IndexError):
                forced[0]
            self.assertEqual(list(forced), [], message)
            with self.assertRaises(StopIteration):
                next(forced._paths)
            self.assertEqual(len(forced), 0, message)
            self.assertFalse(forced, message)
            forced.rewind()
            self.assertEqual(list(forced), [], message)
            self.assertFalse(forced._paths._materialized, message)
        self.check(nw_aligner(epsilon=-1), ("ACGT", "AGT", "+"), (1, 1))

    def test_no_paths_local(self):
        # A local alignment has no paths where no cell scores above 0.
        # Smith-Waterman has none either where it marks no cell as an end
        # point: with an epsilon of -1, a cell must score more than 1 above
        # the highest score before it.
        for aligner, seqA, seqB, score in [
            (sw_aligner(), "AAA", "CCC", 0),
            (sw_aligner(epsilon=-1), "ACGT", "ACGT", 4),
            (gotoh_local_aligner(), "AAA", "CCC", 0),
        ]:
            message = f"{aligner.algorithm} {seqA} {seqB} {aligner.epsilon}"
            expected = aligner.align(seqA, seqB)
            forced = forced_align(aligner, seqA, seqB, budgets=(1, 1))
            self.assertIs(type(forced._paths), LinearPaths, message)
            self.assertEqual(expected.score, score, message)
            self.assertEqual(forced.score, score, message)
            self.assertEqual(len(expected), 0, message)
            self.assertEqual(list(expected), [], message)
            with self.assertRaises(IndexError):
                forced[0]
            self.assertEqual(list(forced), [], message)
            self.assertEqual(len(forced), 0, message)
            self.assertFalse(forced, message)
            self.assertFalse(forced._paths._materialized, message)

    def test_unreached(self):
        # GOTOH_LOCAL_ALIGN's clean-up removes the traces to matrices that no
        # start point reaches, and any end point there.  An infinite gap
        # extension score made every Ix infinite from row 0, where no path
        # can start, and so the only end point, the first infinite M, was
        # reached only from there; with NaN scores as well, the walk met an
        # Ix or Iy without a trace inside the block.  The walk then found
        # the trace it took removed, and asked the full matrix, which had no
        # paths.  Local mode now rejects these positive gap scores, in
        # linear space as well.  A fuzz of Gotoh local aligners (100,000 draws),
        # with match and mismatch scores in {1, 0, -1, inf, -inf, NaN}, gap
        # scores in {-2, -1, 0, -inf, NaN} and epsilons including 0, ±inf
        # and NaN, found no other way into that walk; it stays, for safety.
        unreached = gotoh_local_aligner(open_gap_score=-1, extend_gap_score=math.inf)
        nan_Ix = gotoh_local_aligner(
            match_score=math.nan,
            open_deletion_score=-1,
            extend_deletion_score=math.inf,
        )
        nan_Iy = gotoh_local_aligner(
            mismatch_score=math.nan,
            open_insertion_score=-1,
            extend_insertion_score=math.inf,
            open_deletion_score=-math.inf,
            extend_deletion_score=math.nan,
        )
        for aligner, seqA, seqB, name in (
            (unreached, "CC", "CA", "extend_internal_insertion_score"),
            (
                unreached,
                "AGAACTTATCGCTTGACGTAAC",
                "ACGTTAAC",
                "extend_internal_insertion_score",
            ),
            (nan_Ix, "ACCCC", "AA", "extend_internal_deletion_score"),
            (nan_Iy, "CA", "ACCCA", "extend_internal_insertion_score"),
        ):
            message = f"^gap scores must be 0 or less in local mode, but {name} is inf$"
            with self.assertRaisesRegex(ValueError, message):
                aligner.align(seqA, seqB)
            with self.assertRaisesRegex(ValueError, message):
                forced_align(aligner, seqA, seqB, budgets=(1, 1))
            with self.assertRaisesRegex(ValueError, message):
                aligner.score(seqA, seqB)

    def test_gotoh_local_score_epsilon(self):
        # GOTOH_LOCAL_ALIGN raised its highest score only by more than
        # epsilon, so it stayed 0 while cells scoring epsilon were end
        # points.  Its PathGenerator used to start exhausted then, so that
        # [0] raised IndexError while len() and iteration found the paths,
        # also in Biopython 1.88.  Now all three agree, in both tracebacks,
        # and the score is theirs, as score() has it.
        aligner = gotoh_local_aligner(epsilon=1)
        for seqA, seqB, score, paths in (
            ("A", "A", 1, [[[0, 1], [0, 1]]]),
            ("GAT", "TA", 1, [[[1, 2], [1, 2]], [[2, 3], [0, 1]]]),
            ("A", "C", 0, []),
        ):
            message = f"{seqA} {seqB}"
            expected = aligner.align(seqA, seqB)
            forced = forced_align(aligner, seqA, seqB, budgets=(1, 1))
            self.assertIs(type(forced._paths), LinearPaths, message)
            self.assertEqual(aligner.score(seqA, seqB), score, message)
            for alignments in (expected, forced):
                self.assertEqual(alignments.score, score, message)
                self.assertEqual(len(alignments), len(paths), message)
                for index, path in enumerate(paths):
                    self.assertEqual(outcome(alignments, index), path, message)
                with self.assertRaises(IndexError):
                    alignments[len(paths)]
                found = [alignment.coordinates.tolist() for alignment in alignments]
                self.assertEqual(found, paths, message)
            # [0] before anything else, then len(), in linear space.
            forced = forced_align(aligner, seqA, seqB, budgets=(1, 1))
            if paths:
                self.assertEqual(outcome(forced, 0), paths[0], message)
            self.assertEqual(len(forced), len(paths), message)
            self.assertFalse(forced._paths._materialized, message)

    def test_sw_upstream(self):
        # Upstream defects of the full matrix, also in Biopython 1.86 and
        # 1.88, that the linear-space traceback copied.  With a positive gap
        # score the highest score could be that of a path ending in a gap,
        # which no alignment ends in, and an end point found before it rose
        # stayed one, so that the first path of "AAA" and "CAACC" scored 1,
        # and the PathGenerator then failed on an end point left without a
        # diagonal trace.  Local mode now rejects a positive gap score, in
        # linear space as well.
        aligner = sw_aligner(gap_score=0.5)
        message = "gap scores must be 0 or less in local mode, but "
        message += "open_internal_insertion_score is 0.5"
        for seqA, seqB in (("AAA", "CAACC"), ("ACC", "CACAA")):
            with self.assertRaisesRegex(ValueError, f"^{message}$"):
                aligner.align(seqA, seqB)
            with self.assertRaisesRegex(ValueError, f"^{message}$"):
                forced_align(aligner, seqA, seqB, budgets=(1, 1))
            with self.assertRaisesRegex(ValueError, f"^{message}$"):
                aligner.score(seqA, seqB)
        # With an epsilon of the match score, the highest score rises by no
        # more than epsilon at a time, and the end points found on the way
        # stayed end points: the paths each scored 1 here (three of them
        # for "CC"), and blocked the path scoring 2.
        aligner = sw_aligner(epsilon=1)
        for seqA, seqB in (("CC", "CC"), ("AC", "AC")):
            message = f"{seqA} {seqB}"
            expected = aligner.align(seqA, seqB)
            forced = forced_align(aligner, seqA, seqB, budgets=(1, 1))
            for alignments in (expected, forced):
                self.assertEqual(alignments.score, 2, message)
                self.assertEqual(len(alignments), 1, message)
                self.assertEqual(
                    coordinates(alignments, 2), [[[0, 2], [0, 2]]], message
                )
        # The number of paths overflowed, but adding an overflowed count to
        # a count that had not gave a wrong count: 13274987183379644 here, on
        # Biopython 1.88.
        aligner = sw_aligner(gap_score=0)
        seqA = (
            "AACAAACAAAACCACCCAACCCACACCCCCCAACCAAACCACAACACAACACAACCACCAACAAAC"
            "CCCCCACCCCAAAACAAACCAAACACCCCCACCCAACCCCCACCCCACACACCCACACAAAACCCC"
            "CAACCAACCCCCCCCCAACCCC"
        )
        seqB = (
            "AACAACCAACAAAACACAACCCCCACACACAACCCCAAAAAAACCCAAAACCCAACACACACACAA"
            "CCACACACCCCCCAAAAAACACCACAACACCACCACCACACCAAACAACACAACCCAAAACCAAAC"
            "ACACCAAACAAAAAACACAACA"
        )
        expected = aligner.align(seqA, seqB)
        forced = forced_align(aligner, seqA, seqB, budgets=(1, 1))
        overflow = f"number of optimal alignments is larger than {sys.maxsize}"
        self.assertEqual(count(expected), overflow)
        self.assertEqual(count(forced), overflow)
        self.assertEqual(forced.score, expected.score)
        paths = coordinates(expected, 100)
        self.assertEqual(len(set(map(str, paths))), 100)
        self.assertEqual(coordinates(forced, 100), paths)

    def test_interior_start(self):
        # A path can start at a start point inside the matrix, as the second
        # does here, after the mismatch of A with B.  Both counts seed such a
        # start point with one path, so both count the two optimal local
        # alignments, as enumerating every local alignment gives.
        paths = [[[0, 1], [1, 2]], [[1, 2], [1, 2]]]
        for aligner in (sw_aligner(), gotoh_local_aligner()):
            message = aligner.algorithm
            expected = aligner.align("AA", "BA")
            forced = forced_align(aligner, "AA", "BA", budgets=(1, 1))
            for alignments in (expected, forced):
                self.assertEqual(alignments.score, 1, message)
                self.assertEqual(len(alignments), 2, message)
            self.assertFalse(forced._paths._materialized, message)
            for alignments in (expected, forced):
                self.assertEqual(coordinates(alignments, 3), paths, message)

    def test_overflow(self):
        # TestOverflowError of test_pairwise_aligner, in linear space; the
        # same with affine gap scores; and a count just short of it.
        seqA = SeqIO.read(support.DATA / "Align" / "bsubtilis.fa", "fasta").seq
        seqB = SeqIO.read(support.DATA / "Align" / "ecoli.fa", "fasta").seq
        for aligner, alignments in (
            (
                Align.PairwiseAligner(gap_score=0),
                f">{sys.maxsize} alignments; score=1286",
            ),
            (
                Align.PairwiseAligner(open_gap_score=-0.5, extend_gap_score=0),
                f">{sys.maxsize} alignments; score=1228.5",
            ),
            (
                Align.PairwiseAligner(open_gap_score=-1, extend_gap_score=0),
                "5650532794368000 alignments; score=1201",
            ),
        ):
            expected = aligner.align(seqA, seqB)
            forced = forced_align(aligner, seqA, seqB)
            self.assertEqual(
                repr(forced),
                f"<PairwiseAlignments object ({alignments}) at {hex(id(forced))}>",
            )
            self.assertEqual(count(forced), count(expected))
            first = expected[0].coordinates
            self.assertTrue(np.array_equal(forced[0].coordinates, first))
            self.assertFalse(forced._paths._materialized)
        # Smith-Waterman, overflowing, and a count just short of it.
        aligner = sw_aligner()
        for n, alignments in (
            (None, f">{sys.maxsize} alignments; score=945"),
            (700, "4314952315699200 alignments; score=380"),
        ):
            expected = aligner.align(seqA[:n], seqB[:n])
            forced = forced_align(aligner, seqA[:n], seqB[:n])
            self.assertEqual(
                repr(forced),
                f"<PairwiseAlignments object ({alignments}) at {hex(id(forced))}>",
            )
            self.assertEqual(count(forced), count(expected))
            first = expected[0].coordinates
            self.assertTrue(np.array_equal(forced[0].coordinates, first))
            self.assertFalse(forced._paths._materialized)
        # Gotoh local, overflowing, and a count just short of it.
        for open_gap_score, alignments in (
            (-0.5, f">{sys.maxsize} alignments; score=1228.5"),
            (-1, "11301065588736000 alignments; score=1201"),
        ):
            aligner = gotoh_local_aligner(
                mismatch_score=0, open_gap_score=open_gap_score, extend_gap_score=0
            )
            expected = aligner.align(seqA, seqB)
            forced = forced_align(aligner, seqA, seqB)
            self.assertEqual(
                repr(forced),
                f"<PairwiseAlignments object ({alignments}) at {hex(id(forced))}>",
            )
            self.assertEqual(count(forced), count(expected))
            first = expected[0].coordinates
            self.assertTrue(np.array_equal(forced[0].coordinates, first))
            self.assertFalse(forced._paths._materialized)
        # Counts after an overflow, which a count that had not overflowed
        # used to absorb (see test_count_past_overflow in
        # test_pairwise_aligner): 26 by 26 gives the largest central
        # Delannoy number below sys.maxsize, and 37 by 37 a wrong count.
        aligner = nw_aligner(match_score=0, mismatch_score=0, gap_score=0)
        overflow = f"number of optimal alignments is larger than {sys.maxsize}"
        for n, alignments in (
            (26, 8970232353223635949),
            (27, overflow),
            (37, overflow),
        ):
            forced = forced_align(aligner, "A" * n, "C" * n)
            self.assertEqual(count(forced), alignments, n)
            self.assertEqual(count(aligner.align("A" * n, "C" * n)), alignments, n)
            self.assertFalse(forced._paths._materialized)

    @unittest.skipUnless(SLOW, "set BIOPAITHON_SLOW_TESTS=1 to run")
    def test_soak(self):
        self.run_cases(seed=4589, count=50000, lengths=(1, 60))
        self.run_cases(seed=3407, count=300, lengths=(100, 1500))
        self.run_cases(seed=4590, count=50000, lengths=(1, 60), algorithm=GOTOH)
        self.run_cases(seed=3408, count=300, lengths=(100, 1500), algorithm=GOTOH)
        self.run_cases(seed=4591, count=50000, lengths=(1, 60), algorithm=SW)
        self.run_cases(seed=3409, count=300, lengths=(100, 1500), algorithm=SW)
        self.run_cases(seed=4592, count=50000, lengths=(1, 60), algorithm=GOTOH_LOCAL)
        self.run_cases(seed=3410, count=300, lengths=(100, 1500), algorithm=GOTOH_LOCAL)

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
        for aligner in (
            nw_aligner(),
            gotoh_aligner(),
            sw_aligner(),
            gotoh_local_aligner(),
        ):
            alignments = forced_align(aligner, seqA, seqB, budgets=DEFAULT[1:])
            alignment = alignments[0]
            self.assertEqual(alignment.counts(aligner).score, alignments.score)
            self.assertEqual(alignments.score, aligner.score(seqA, seqB))
            if aligner.mode == "global":
                end = alignment.coordinates[:, -1].tolist()
                self.assertEqual(end, [len(seqA), len(seqB)])

    @unittest.skipUnless(SLOW, "set BIOPAITHON_SLOW_TESTS=1 to run")
    def test_read_against_genome(self):
        # A 1,000-letter read against 2.2 million letters, with the default
        # limits: over INT_MAX cells, and a 2.2 GB full matrix (4.4 GB for
        # Gotoh).  Aligned locally, the read is found where it came from.
        rng = random.Random(2200000)
        seqA = random_dna(rng, 2200000)
        seqB = list(seqA[1100000:1101000])
        for _ in range(20):
            seqB[rng.randrange(len(seqB))] = rng.choice("ACGT")
        seqB = "".join(seqB)
        self.assertGreater((len(seqA) + 1) * (len(seqB) + 1), 2**31 - 1)
        for aligner in (
            nw_aligner(),
            gotoh_aligner(),
            sw_aligner(),
            gotoh_local_aligner(),
        ):
            alignments = aligner.align(seqA, seqB)
            self.assertIs(type(alignments._paths), LinearPaths)
            alignment = alignments[0]
            self.assertEqual(alignment.counts(aligner).score, alignments.score)
            self.assertEqual(alignments.score, aligner.score(seqA, seqB))
            self.assertFalse(alignments._paths._materialized)
            if aligner.mode == "local":
                start, end = alignment.coordinates[0, [0, -1]].tolist()
                self.assertLessEqual(1100000, start)
                self.assertLessEqual(end, 1101000)
                self.assertGreater(end - start, 900)


class TestWatermanSmithBeyerLocal(unittest.TestCase):
    """The Waterman-Smith-Beyer local full matrix against brute force.

    It has no linear-space traceback, but shares the brute force here.
    """

    def test_brute_force(self):
        # As test_gotoh_local_brute_force, with gap score functions of the
        # length of the gap, affine or not.
        rng = random.Random(2032)
        values = [-2, -1, -0.5, 0, 0.5, 1]
        # local mode rejects positive gap scores
        gap_values = [value for value in values if value <= 0]
        for _ in range(1500):
            alphabet = rng.choice(["AC", "ACG"])
            seqA = "".join(rng.choice(alphabet) for _ in range(rng.randint(1, 5)))
            seqB = "".join(rng.choice(alphabet) for _ in range(rng.randint(1, 5)))
            if rng.random() < 0.5:
                match = rng.choice([1, 2, 0.5])
                mismatch = rng.choice([-1, 0, -0.5])
                scores = {"match_score": match, "mismatch_score": mismatch}

                def substitution(a, b, match=match, mismatch=mismatch):
                    return match if a == b else mismatch

            else:
                matrix = substitution_matrices.Array(alphabet, dims=2)
                matrix[:, :] = [
                    [rng.choice(values) for _ in alphabet] for _ in alphabet
                ]
                scores = {"substitution_matrix": matrix}

                def substitution(a, b, matrix=matrix):
                    return matrix[a, b]

            functions = []
            for _ in range(2):
                if rng.random() < 0.5:
                    opening, extension = rng.choice(gap_values), rng.choice(gap_values)
                    table = [opening + extension * k for k in range(5)]
                else:
                    table = [rng.choice(gap_values) for _ in range(5)]
                functions.append(table)
            if rng.random() < 0.7:
                functions[1] = functions[0]
            insertion_table, deletion_table = functions

            def insertion(n, table=insertion_table):
                return table[n - 1]

            def deletion(n, table=deletion_table):
                return table[n - 1]

            aligner = Align.PairwiseAligner(
                mode="local",
                insertion_score=lambda i, n: insertion(n),
                deletion_score=lambda i, n: deletion(n),
                **scores,
            )
            self.assertEqual(aligner.algorithm, WSB_LOCAL)
            # The scores are dyadic, so they add up exactly.
            aligner.epsilon = rng.choice([1e-6, 0.25, 0.5, 1, 2, 5])
            message = f"{aligner} {seqA} {seqB} {insertion_table} {deletion_table}"
            best, paths = brute_force_wsb_local(
                seqA, seqB, substitution, insertion, deletion
            )
            self.assertEqual(aligner.score(seqA, seqB), best, message)
            tied = tied_wsb_local(
                seqA, seqB, substitution, insertion, deletion, aligner.epsilon
            )
            if aligner.epsilon < 0.5:
                self.assertEqual(tied, (best, paths), message)
            best, paths = tied
            alignments = aligner.align(seqA, seqB)
            self.assertEqual(alignments.score, best, message)
            self.assertEqual(len(alignments), len(paths), message)
            found = [
                tuple(map(tuple, coordinates))
                for coordinates in coordinates(alignments, len(paths) + 1)
            ]
            self.assertEqual(len(found), len(paths), message)
            self.assertEqual(set(found), paths, message)
            indexed = [
                tuple(map(tuple, alignments[k].coordinates.tolist()))
                for k in range(len(paths))
            ]
            self.assertEqual(indexed, found, message)


class TestSnapshot(unittest.TestCase):
    """Changes made after align() do not change the alignments."""

    def check(self, aligner, seqA, seqB, change):
        expected = aligner.align(seqA, seqB)
        reference = (expected[0].coordinates, expected[1].coordinates, len(expected))
        self.assertGreater(reference[2], 2)
        forced = forced_align(aligner, seqA, seqB, budgets=(1, 1))
        self.assertIs(type(forced._paths), LinearPaths, aligner.algorithm)
        change()
        self.assertTrue(np.array_equal(forced[0].coordinates, reference[0]))
        self.assertTrue(np.array_equal(forced[1].coordinates, reference[1]))
        self.assertEqual(len(forced), reference[2])

    def test_scores(self):
        seqA = "GAACTNGACGTTANCC"
        seqB = "GACNTGACGTAACCA"

        def change(aligner):
            aligner.match_score = -3
            aligner.mismatch_score = 2
            aligner.gap_score = 0.5
            aligner.open_internal_deletion_score = -1.5
            aligner.epsilon = 10
            aligner.wildcard = "A"

        for aligner in (
            nw_aligner(wildcard="N"),
            gotoh_aligner(
                wildcard="N", mismatch_score=0, open_gap_score=-0.5, extend_gap_score=0
            ),
            sw_aligner(wildcard="N", mismatch_score=0, gap_score=0),
            gotoh_local_aligner(
                wildcard="N", mismatch_score=0, open_gap_score=-0.5, extend_gap_score=0
            ),
        ):
            self.check(aligner, seqA, seqB, functools.partial(change, aligner))

    def test_array(self):
        seqA = "GAACTTGACGTTAGCC"
        seqB = "GACATGACGTAACCA"

        def change(matrix):
            for letter in "ACGT":
                matrix[letter, letter] = -10

        for make in (nw_aligner, gotoh_aligner, sw_aligner, gotoh_local_aligner):
            matrix = substitution_matrices.load("NUC.4.4")
            aligner = make(substitution_matrix=matrix)
            self.check(aligner, seqA, seqB, functools.partial(change, matrix))

    def test_caller_arrays(self):
        def change(seqA, seqB, matrix):
            seqA[:] = 0
            seqB[:] = 1
            matrix[:, :] = 5

        gotoh_local = functools.partial(
            gotoh_local_aligner, open_gap_score=-1, extend_gap_score=-0.5
        )
        for make in (nw_aligner, gotoh_aligner, sw_aligner, gotoh_local):
            matrix = 2 * np.identity(4) - 1
            aligner = make(substitution_matrix=matrix)
            seqA = np.array(["ACGT".index(c) for c in "GAACTTGACGTTAGCCTA"], np.int32)
            seqB = np.array(["ACGT".index(c) for c in "GACATGACGTAACCA"], np.int32)
            change_all = functools.partial(change, seqA, seqB, matrix)
            self.check(aligner, seqA, seqB, change_all)


class TestThreads(unittest.TestCase):
    """Concurrent next(), len() and reset() on one paths object."""

    def test_threads(self):
        for aligner in (
            nw_aligner(mismatch_score=0, gap_score=0),
            gotoh_aligner(open_gap_score=-0.5, extend_gap_score=-1),
            sw_aligner(mismatch_score=0, gap_score=0),
            gotoh_local_aligner(open_gap_score=-0.5, extend_gap_score=-1),
        ):
            with self.subTest(algorithm=aligner.algorithm):
                self.check(aligner)

    def check(self, aligner):
        seqA = "GAACTTGACGTTAGCCTA"
        seqB = "GACATGACGTAACCA"
        expected = aligner.align(seqA, seqB)
        length = len(expected)
        self.assertGreater(length, 10)
        paths = set(expected._paths)
        self.assertEqual(len(paths), length)
        shared = forced_align(aligner, seqA, seqB, budgets=(1, 1))._paths
        self.assertIs(type(shared), LinearPaths)
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
        cls.gotoh_aligner = gotoh_aligner()
        cls.gotoh_expected = cls.gotoh_aligner.align(*cls.seqs)
        cls.sw_aligner = sw_aligner()
        cls.sw_expected = cls.sw_aligner.align(*cls.seqs)
        # The best local alignment of two random sequences is short, so its
        # first path took about 5 ms, close enough to the three 1 ms timer
        # ticks that it could finish first on a fast machine.  A query that
        # matches each copy of a tandem repeat gives two long alignments.
        block = random_dna(rng, 2000)
        mutated = "".join(
            rng.choice("ACGT") if rng.random() < 0.15 else c for c in block
        )
        cls.gotoh_local_seqs = (block + block, mutated)
        cls.gotoh_local_aligner = gotoh_local_aligner()
        cls.gotoh_local_expected = cls.gotoh_local_aligner.align(*cls.gotoh_local_seqs)

    @contextlib.contextmanager
    def alarms(self, handler, delay, interval=0):
        """Run handler on SIGALRM from a timer for the duration of a with block."""
        previous = signal.signal(signal.SIGALRM, handler)
        signal.setitimer(signal.ITIMER_REAL, delay, interval)
        try:
            yield
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
            # Discard a SIGALRM still pending for another thread, such as a
            # pytest-xdist worker's: under the default action, restored next,
            # it would end the process.
            signal.signal(signal.SIGALRM, signal.SIG_IGN)
            signal.signal(signal.SIGALRM, previous)

    def interrupt(self, paths, returned, call=next):
        """Call call(paths) with a timer raising KeyboardInterrupt."""
        calls = 0

        def handler(signum, frame):
            nonlocal calls
            calls += 1
            if calls == 3:
                raise KeyboardInterrupt

        with self.alarms(handler, 0.001, 0.001):
            returned.append(call(paths))

    def test_interrupt(self):
        stages = (("first path", 5), ("full matrix", 5), ("count", 3))
        self.check_interrupt(self.aligner, self.expected, stages)

    def test_interrupt_gotoh(self):
        stages = (("first path", 5), ("full matrix", 5), ("count", 3))
        self.check_interrupt(self.gotoh_aligner, self.gotoh_expected, stages)

    def test_interrupt_sw(self):
        stages = (("first path", 5), ("full matrix", 5), ("count", 3))
        self.check_interrupt(self.sw_aligner, self.sw_expected, stages)

    def test_interrupt_gotoh_local(self):
        stages = (("first path", 5), ("full matrix", 5), ("count", 3))
        aligner, expected = self.gotoh_local_aligner, self.gotoh_local_expected
        self.check_interrupt(aligner, expected, stages, self.gotoh_local_seqs)

    def check_interrupt(self, aligner, expected, stages, seqs=None):
        if seqs is None:
            seqs = self.seqs
        reference = (expected[0].coordinates, expected[1].coordinates, count(expected))
        for stage, repeats in stages:
            for _ in range(repeats):
                alignments = forced_align(aligner, *seqs, budgets=DEFAULT[1:])
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
        #
        # Python runs the handler only between bytecodes or where the C
        # code checks for signals, so a single tick can land before next()
        # enters C, where len() would count every path; under load a 2 ms
        # tick landed in the re.compile of assertRaisesRegex.  So the
        # timer ticks every millisecond and each tick first probes with
        # reset(), which outside next() is cheap and harmless; the first
        # tick inside next() fails, stops the timer so no tick nests in
        # it, checks that len() fails the same way, and its error ends
        # next().
        def handler(signum, frame):
            try:
                paths.reset()
            except RuntimeError:
                signal.setitimer(signal.ITIMER_REAL, 0)
                with self.assertRaisesRegex(RuntimeError, "being computed"):
                    len(paths)
                raise

        alignments = forced_align(self.aligner, *self.seqs)
        paths = alignments._paths
        with self.alarms(handler, 0.001, 0.001):
            with self.assertRaisesRegex(RuntimeError, "being computed"):
                next(paths)
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
calls = {"tests": 0, NW: 0, GOTOH: 0, SW: 0, GOTOH_LOCAL: 0}


def checked_align(self, *args, **kwargs):
    """Check [0] against the full matrix, then align in linear space."""
    with align_lock:
        forced = original_align(self, *args, **kwargs)
        algorithm = self.algorithm
        if algorithm not in (NW, GOTOH, SW, GOTOH_LOCAL):
            return forced
        # Every Needleman-Wunsch, Smith-Waterman and Gotoh alignment takes
        # the linear-space traceback.
        assert type(forced._paths) is LinearPaths, type(forced._paths)
        calls[algorithm] += 1
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


# TestMatrixAllocationGuards never reaches the linear-space traceback: its
# alignments fail before routing, or are Waterman-Smith-Beyer or FOGSAA.  So a
# forced copy would only run it a second time.
NOT_FORCED = {"TestMatrixAllocationGuards"}


def forced_classes(module):
    """Return a forced re-run class for each test case class of a module."""
    classes = {}
    for name, value in vars(module).items():
        if (
            isinstance(value, type)
            and issubclass(value, unittest.TestCase)
            and value.__module__ == module.__name__
            and name not in NOT_FORCED
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
        self.assertGreater(calls[NW], 50)
        self.assertGreater(calls[GOTOH], 50)
        self.assertGreater(calls[SW], 40)
        self.assertGreater(calls[GOTOH_LOCAL], 50)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
