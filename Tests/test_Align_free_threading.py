# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests of the Bio.Align C extensions called from several threads at once.

These matter only without the GIL, so they run only on a free-threaded
build of Python. Each runs in a child process, so that a crash or a hang
fails the test rather than taking down the whole test run. The child
skips if the GIL is enabled, as it would be if Bio.Align imported an
extension that has not declared it can run without the GIL.
"""

import os
import subprocess
import sys
import sysconfig
import unittest

try:
    import numpy  # noqa: F401
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align."
    ) from None


# Six threads iterate over one path generator while a seventh keeps
# resetting it and taking its length, for each kind of path generator.
SHARED_PATH_GENERATOR = """\
import sys
from concurrent.futures import ThreadPoolExecutor
from threading import Barrier, Event

from Bio.Align import CodonAligner, PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

if sys._is_gil_enabled():
    print("SKIP: the GIL is enabled")
    sys.exit()


def check(aligner, seqA, seqB):
    # Every path a thread gets must be one of the optimal paths, whatever
    # the other threads do to the shared path generator.
    expected = set(aligner.align(seqA, seqB)._paths)
    paths = aligner.align(seqA, seqB)._paths
    barrier = Barrier(7)
    done = Event()

    def iterate():
        barrier.wait()
        for _ in range(3000):
            try:
                path = next(paths)
            except StopIteration:
                continue
            if path not in expected:
                raise AssertionError(f"unexpected path {path}")

    def reset():
        barrier.wait()
        while not done.is_set():
            paths.reset()
            if len(paths) != len(expected):
                raise AssertionError(f"{len(paths)} paths, expected {len(expected)}")

    with ThreadPoolExecutor(max_workers=7) as executor:
        resetter = executor.submit(reset)
        try:
            iterators = [executor.submit(iterate) for _ in range(6)]
            for future in iterators:
                future.result()
        finally:
            done.set()
        resetter.result()


aligner = PairwiseAligner(mismatch_score=-1, gap_score=-1)  # Needleman-Wunsch
check(aligner, "GAACTACTACTACTGA" * 2, "GACTACTGA" * 2)
aligner = PairwiseAligner(mode="local", open_gap_score=-1, extend_gap_score=-0.5)
check(aligner, "A" * 30, "A" * 10)  # Gotoh
aligner = PairwiseAligner(gap_score=lambda start, length: -1 - 0.5 * (length - 1))
check(aligner, "A" * 30, "A" * 10)  # Waterman-Smith-Beyer
aligner = CodonAligner()
aligner.frameshift_score = 0.0
protein = SeqRecord(Seq("FKKKF"), id="protein")
dna = SeqRecord(Seq("TTTAAAAAAAAAATTT"), id="dna")
check(aligner, protein, dna)
print("OK")
"""

# Six threads score and align a shared NumPy-encoded query while two others
# keep replacing the aligner's substitution matrix and gap score, the
# latter through each way of setting an attribute.
SCORE_WHILE_RECONFIGURED = """\
import sys
from concurrent.futures import ThreadPoolExecutor
from threading import Barrier, Event

import numpy as np

from Bio.Align import PairwiseAligner, _pairwisealigner, substitution_matrices

if sys._is_gil_enabled():
    print("SKIP: the GIL is enabled")
    sys.exit()

names = ("BLOSUM62", "BLOSUM80")


def gap_function(start, length):
    return -10.0 - length


def set_gap_score(aligner, use_function):
    # A new function each time, so that replacing it drops the last
    # reference to the old one; the same goes for each loaded matrix.
    if use_function:
        aligner.gap_score = lambda start, length: gap_function(start, length)
    else:
        aligner.gap_score = -7.0


def encode(text):
    return np.frombuffer(text.encode("utf-32-le"), np.int32).copy()


query = encode("MKTAYIAKQRQISFVKSHFSRQLEE")
query_copy = query.copy()
targets = [encode("MKTAYIAKQWRQISFVKSHFRQLEEP"), encode("MKAYIAKQRQISFVHFSRQLE")]

# The scores of each target under every configuration; a call that sees one
# consistent configuration returns one of these.
aligner = PairwiseAligner()
expected = [set(), set()]
for name in names:
    for use_function in (False, True):
        aligner.substitution_matrix = substitution_matrices.load(name)
        set_gap_score(aligner, use_function)
        for scores, target in zip(expected, targets):
            scores.add(aligner.score(query, target))
            scores.add(aligner.align(query, target).score)

barrier = Barrier(8)
done = Event()


def score(seed):
    rng = np.random.default_rng(seed)
    barrier.wait()
    for _ in range(150):
        k = int(rng.integers(2))
        if rng.random() < 0.5:
            value = aligner.score(query, targets[k])
        else:
            value = aligner.align(query, targets[k]).score
        if value not in expected[k]:
            raise AssertionError(f"score {value} is not one of {expected[k]}")
        str(aligner)


def switch_matrix():
    barrier.wait()
    i = 0
    while not done.is_set():
        i += 1
        aligner.substitution_matrix = substitution_matrices.load(names[i % 2])


descriptor = _pairwisealigner.PairwiseAligner.gap_score


def switch_gap_score():
    barrier.wait()
    i = 0
    while not done.is_set():
        i += 1
        if i % 2:
            value = lambda start, length: gap_function(start, length)
        else:
            value = -7.0
        route = i % 6 // 2
        if route == 0:
            aligner.gap_score = value
        elif route == 1:
            object.__setattr__(aligner, "gap_score", value)
        else:
            descriptor.__set__(aligner, value)
        descriptor.__get__(aligner)


with ThreadPoolExecutor(max_workers=8) as executor:
    switchers = [executor.submit(switch_matrix), executor.submit(switch_gap_score)]
    try:
        scorers = [executor.submit(score, seed) for seed in range(6)]
        for future in scorers:
            future.result()
    finally:
        done.set()
    for future in switchers:
        future.result()
if not np.array_equal(query, query_copy):
    raise AssertionError("the shared query was changed")
print("OK")
"""

# Six threads score and align one shared NumPy-encoded query against the
# same targets with one shared aligner, and must get the serial results.
SCORE_SHARED_QUERY = """\
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from Bio.Align import PairwiseAligner

if sys._is_gil_enabled():
    print("SKIP: the GIL is enabled")
    sys.exit()

rng = np.random.default_rng(0)
letters = np.array(list("ACDEFGHIKLMNPQRSTVWY"))


def random_protein(length):
    return "".join(rng.choice(letters, length))


query = np.frombuffer(random_protein(300).encode("utf-32-le"), np.int32).copy()
query_copy = query.copy()
targets = [random_protein(int(rng.integers(200, 400))) for _ in range(24)]
aligner = PairwiseAligner(scoring="blastp")
expected = [
    (aligner.score(query, target), aligner.align(query, target)[0].coordinates)
    for target in targets
]


def work(seed):
    for k in np.random.default_rng(seed).permutation(len(targets)):
        score = aligner.score(query, targets[k])
        alignment = aligner.align(query, targets[k])[0]
        if score != expected[k][0] or alignment.score != score:
            raise AssertionError(f"target {k}: score {score}, expected {expected[k][0]}")
        if not np.array_equal(alignment.coordinates, expected[k][1]):
            raise AssertionError(f"target {k}: the alignment differs")


with ThreadPoolExecutor(max_workers=6) as executor:
    for future in [executor.submit(work, seed) for seed in range(12)]:
        future.result()
if not np.array_equal(query, query_copy):
    raise AssertionError("the shared query was changed")
print("OK")
"""

# Four threads feed lines to one printed-alignment parser while two others
# keep reading its shape.
SHARED_PARSER = """\
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from threading import Barrier, Event

import numpy as np

from Bio.Align._aligncore import PrintedAlignmentParser

if sys._is_gil_enabled():
    print("SKIP: the GIL is enabled")
    sys.exit()

rng = np.random.default_rng(0)
nlines = 400
chunks = [
    [bytes(rng.choice(list(b"ACGT---"), 80).tolist()) for _ in range(nlines)]
    for _ in range(4)
]
parser = PrintedAlignmentParser()
barrier = Barrier(6)
done = Event()


def feed(lines):
    barrier.wait()
    for line in lines:
        nbytes, sequence = parser.feed(line)
        if nbytes != len(line) or sequence != line.replace(b"-", b""):
            raise AssertionError(f"feed({line}) returned {nbytes}, {sequence}")


def shape():
    barrier.wait()
    previous = 0
    while not done.is_set():
        n, k = parser.shape
        if not previous <= n <= 4 * nlines or not 1 <= k <= 81:
            raise AssertionError(f"shape ({n}, {k}) after {previous} lines")
        previous = n


with ThreadPoolExecutor(max_workers=6) as executor:
    readers = [executor.submit(shape) for _ in range(2)]
    try:
        feeders = [executor.submit(feed, lines) for lines in chunks]
        for future in feeders:
            future.result()
    finally:
        done.set()
    for future in readers:
        future.result()

# Every line is in the parser exactly once, in some order.
n, k = parser.shape
if n != 4 * nlines:
    raise AssertionError(f"{n} lines in the parser")
coordinates = np.zeros((n, k), np.intp)
parser.fill(coordinates)
fed = Counter(len(line.replace(b"-", b"")) for lines in chunks for line in lines)
if Counter(coordinates[:, -1].tolist()) != fed:
    raise AssertionError("lines were lost or changed")
print("OK")
"""

# Six threads count an alignment with one aligner while two others keep
# replacing the aligner's substitution matrix and gap score.
COUNTS_WHILE_RECONFIGURED = """\
import sys
from concurrent.futures import ThreadPoolExecutor
from threading import Barrier, Event

from Bio.Align import PairwiseAligner, substitution_matrices

if sys._is_gil_enabled():
    print("SKIP: the GIL is enabled")
    sys.exit()

names = ("BLOSUM62", "BLOSUM80")


def gap_function(start, length):
    return -10.0 - length


def set_gap_score(aligner, use_function):
    # A new function each time, so that replacing it drops the last
    # reference to the old one; the same goes for each loaded matrix.
    if use_function:
        aligner.gap_score = lambda start, length: gap_function(start, length)
    else:
        aligner.gap_score = -7.0


aligner = PairwiseAligner(scoring="blastp")
alignment = aligner.align("MKTAYIAKQRQISFVKSHFSRQLEE", "MKTAYIAKQWRQISFVKSHFRQLEEP")[0]

# The counts under every configuration; a call that sees one consistent
# configuration returns one of these.
expected = set()
for name in names:
    for use_function in (False, True):
        aligner.substitution_matrix = substitution_matrices.load(name)
        set_gap_score(aligner, use_function)
        counts = alignment.counts(aligner)
        expected.add((counts.substitution_score, counts.gap_score))

barrier = Barrier(8)
done = Event()


def count():
    barrier.wait()
    for _ in range(300):
        counts = alignment.counts(aligner)
        observed = (counts.substitution_score, counts.gap_score)
        if observed not in expected:
            raise AssertionError(f"counts {observed} are not one of {expected}")


def switch_matrix():
    barrier.wait()
    i = 0
    while not done.is_set():
        i += 1
        aligner.substitution_matrix = substitution_matrices.load(names[i % 2])


def switch_gap_score():
    barrier.wait()
    i = 0
    while not done.is_set():
        i += 1
        set_gap_score(aligner, i % 2)


with ThreadPoolExecutor(max_workers=8) as executor:
    switchers = [executor.submit(switch_matrix), executor.submit(switch_gap_score)]
    try:
        counters = [executor.submit(count) for _ in range(6)]
        for future in counters:
            future.result()
    finally:
        done.set()
    for future in switchers:
        future.result()
print("OK")
"""

# Four threads count an alignment of lazy sequences, through the private
# AlignmentCounts constructor, while two others keep replacing the
# sequences in the list they all share and changing its length.
COUNTS_SHARED_LIST = """\
import sys
from concurrent.futures import ThreadPoolExecutor
from threading import Barrier, Event

import numpy as np

from Bio.Align import _alignmentcounts
from Bio.Seq import SequenceDataAbstractBaseClass

if sys._is_gil_enabled():
    print("SKIP: the GIL is enabled")
    sys.exit()


class LazyData(SequenceDataAbstractBaseClass):
    # Sequence data read on demand, through the sequence protocol.

    def __init__(self, data):
        self.data = data
        super().__init__()

    def __len__(self):
        return len(self.data)

    def __getitem__(self, key):
        return self.data[key]


target = b"ACGTACGTAC" * 10
query = b"ACGTTCGTAC" * 10
sequences = [LazyData(target), LazyData(query)]
coordinates = np.array([np.arange(0, 101, 5), np.arange(0, 101, 5)])
strands = np.zeros(2, bool)
barrier = Barrier(6)
done = Event()


def count():
    barrier.wait()
    successes = 0
    for _ in range(2000):
        try:
            counts = _alignmentcounts.AlignmentCounts(sequences, coordinates, strands)
        except ValueError as exception:
            # The list had more than two sequences when copied.
            if str(exception) != (
                "number of rows in coordinates must equal the number of sequences"
            ):
                raise
            continue
        if (counts.identities, counts.mismatches) != (90, 10):
            raise AssertionError(f"{counts.identities} identities")
        successes += 1
    # About half the calls succeed; none would mean the counts went unchecked.
    if successes == 0:
        raise AssertionError("no count succeeded")


def change():
    # Replace the sequences with new objects, freeing the old ones, and
    # change the length of the list.
    barrier.wait()
    while not done.is_set():
        sequences[0] = LazyData(target)
        sequences[1] = LazyData(query)
        sequences.append(LazyData(query))
        del sequences[2]


with ThreadPoolExecutor(max_workers=6) as executor:
    changers = [executor.submit(change) for _ in range(2)]
    try:
        counters = [executor.submit(count) for _ in range(4)]
        for future in counters:
            future.result()
    finally:
        done.set()
    for future in changers:
        future.result()
print("OK")
"""

# Four threads race to set the alphabet of a new substitution matrix.
SET_ALPHABET = """\
import sys
from threading import Barrier, Thread

import numpy as np

from Bio.Align import PairwiseAligner
from Bio.Align.substitution_matrices import Array

if sys._is_gil_enabled():
    print("SKIP: the GIL is enabled")
    sys.exit()

# Alphabets of 1, 2 and 4 bytes per letter, which have mappings of
# different sizes, and a tuple, which has none.
alphabets = (
    "ACGT",
    "".join(map(chr, range(0x391, 0x395))),
    "".join(map(chr, range(0x1D538, 0x1D53C))),
    ("A", "C", "G", "T"),
)
data = np.arange(16.0).reshape(4, 4)
aligner = PairwiseAligner(gap_score=-100)
barrier = Barrier(len(alphabets))

for _ in range(200):
    matrix = data.copy().view(Array)  # no alphabet yet
    won = []

    def set_alphabet(alphabet):
        barrier.wait()
        try:
            matrix.alphabet = alphabet
        except ValueError as exception:
            if str(exception) != "the alphabet has already been set.":
                raise
        else:
            won.append(alphabet)

    threads = [Thread(target=set_alphabet, args=(alphabet,)) for alphabet in alphabets]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    if len(won) != 1:
        raise AssertionError(f"{len(won)} threads set the alphabet")
    alphabet = won[0]
    if matrix.alphabet is not alphabet:
        raise AssertionError(f"alphabet {matrix.alphabet}, expected {alphabet}")
    # The mapping must belong to the alphabet that was set.
    aligner.substitution_matrix = matrix
    for i, letter1 in enumerate(alphabet):
        for j, letter2 in enumerate(alphabet):
            if isinstance(alphabet, str):  # mapped by the matrix's mapping
                score = aligner.score(letter1, letter2)
            else:  # mapped to indices in Python
                score = aligner.score([letter1], [letter2])
            if score != data[i, j]:
                raise AssertionError(f"{alphabet}: {letter1}, {letter2} scored {score}")
print("OK")
"""

# Four threads keep giving a shared aligner a substitution matrix whose
# __release_buffer__ method sleeps, and then dropping it again by setting
# match_score, mismatch_score or substitution_matrix. Sleeping detaches the
# thread, and so lets another thread take the aligner's lock while a setter
# is releasing the old matrix.
RELEASE_MATRIX_WHILE_RECONFIGURED = """\
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from threading import Barrier

import numpy as np

from Bio.Align import PairwiseAligner

if sys._is_gil_enabled():
    print("SKIP: the GIL is enabled")
    sys.exit()


class Matrix(np.ndarray):
    def __release_buffer__(self, view):
        time.sleep(1e-5)


aligner = PairwiseAligner()
sequence = np.arange(4, dtype=np.int32)
settings = (
    ("match_score", 1.0),
    ("mismatch_score", 0.0),
    ("substitution_matrix", None),
    ("substitution_matrix", np.eye(4)),
)
barrier = Barrier(4)


def change(i):
    barrier.wait()
    for j in range(2000):
        aligner.substitution_matrix = np.eye(4).view(Matrix)
        name, value = settings[(i + j) % len(settings)]
        setattr(aligner, name, value)
        # Every setting scores four matches as 4.0.
        score = aligner.score(sequence, sequence)
        if score != 4.0:
            raise AssertionError(f"score {score} after setting {name}")


with ThreadPoolExecutor(max_workers=4) as executor:
    for future in [executor.submit(change, i) for i in range(4)]:
        future.result()
print("OK")
"""

# Eight threads share the alignments of one linear-space traceback, mostly
# taking their length, which takes and releases their lock each time. A
# thread that held the lock before must not take itself for the holder
# while another thread is taking it.
SHARED_LINEAR_PATHS = """\
import sys
from concurrent.futures import ThreadPoolExecutor
from threading import Barrier

from Bio.Align import PairwiseAligner, _pairwisealigner

if sys._is_gil_enabled():
    print("SKIP: the GIL is enabled")
    sys.exit()

aligner = PairwiseAligner(mismatch_score=0, gap_score=0)  # Needleman-Wunsch
seqA = "GAACTTGACGTTAGCCTA"
seqB = "GACATGACGTAACCA"
expected = set(aligner.align(seqA, seqB)._paths)
previous = _pairwisealigner._set_traceback_limits(-1, 1, 1)  # always linear
try:
    paths = aligner.align(seqA, seqB)._paths
finally:
    _pairwisealigner._set_traceback_limits(*previous)
if type(paths).__name__ != "LinearPaths":
    raise AssertionError(f"got {type(paths).__name__}")
barrier = Barrier(8)


def work(i):
    barrier.wait()
    for j in range(20000):
        if j % 100 == i:
            paths.reset()
        elif j % 10 == 0:
            try:
                path = next(paths)
            except StopIteration:
                continue
            if path not in expected:
                raise AssertionError(f"unexpected path {path}")
        elif len(paths) != len(expected):
            raise AssertionError(f"{len(paths)} paths, expected {len(expected)}")


with ThreadPoolExecutor(max_workers=8) as executor:
    for future in [executor.submit(work, i) for i in range(8)]:
        future.result()
print("OK")
"""


@unittest.skipUnless(
    sysconfig.get_config_var("Py_GIL_DISABLED"), "requires a free-threaded build"
)
class ThreadTests(unittest.TestCase):
    """Bio.Align's C extensions, called from several threads at once."""

    def run_child(self, code, env=None):
        result = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            env=env,
            text=True,
            timeout=300,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        if result.stdout.startswith("SKIP"):
            self.skipTest(result.stdout.strip())
        self.assertEqual(result.stdout.strip(), "OK", result.stderr)

    def test_shared_path_generator(self):
        """Check threads can iterate over, reset and count one set of alignments.

        Each thread must get only optimal alignments.
        """
        self.run_child(SHARED_PATH_GENERATOR)

    def test_score_while_reconfigured(self):
        """Check threads can score while others replace the matrix and gap score.

        Each score must be that of one of the configurations the aligner had,
        whether the gap score is set as an attribute, by object.__setattr__,
        or through its descriptor.
        """
        self.run_child(SCORE_WHILE_RECONFIGURED)

    def test_score_shared_query(self):
        """Check threads scoring one query with one aligner get serial results."""
        self.run_child(SCORE_SHARED_QUERY)

    def test_shared_parser(self):
        """Check threads can feed one printed-alignment parser and read its shape.

        Every line fed must end up in the parser exactly once.
        """
        self.run_child(SHARED_PARSER)

    def test_counts_while_reconfigured(self):
        """Check threads can count while others replace the matrix and gap score.

        Each count must be that of one of the configurations the aligner had.
        """
        self.run_child(COUNTS_WHILE_RECONFIGURED)

    def test_counts_shared_list(self):
        """Check threads can count while others change the list of sequences.

        Each count must be right, or fail cleanly if the list was the wrong
        length when it was read.
        """
        self.run_child(COUNTS_SHARED_LIST)

    def test_set_alphabet(self):
        """Check only one of several threads can set a matrix's alphabet.

        The mapping used for scoring must belong to the alphabet that was set.
        """
        self.run_child(SET_ALPHABET)

    def test_release_matrix_while_reconfigured(self):
        """Check threads can drop a matrix whose buffer release lets others in.

        Each export of the matrix must be released once. PYTHONMALLOC=debug
        overwrites freed memory, so that releasing one twice crashes.
        """
        env = dict(os.environ, PYTHONMALLOC="debug")
        self.run_child(RELEASE_MATRIX_WHILE_RECONFIGURED, env=env)

    def test_shared_linear_paths(self):
        """Check threads can count, iterate over and reset linear-space alignments.

        No thread may take itself for the one computing them.
        """
        self.run_child(SHARED_LINEAR_PATHS)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
