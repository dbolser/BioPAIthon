# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""What a downstream ``mypy --strict`` run sees of the Bio.Align C extension stubs.

PairwiseAligner inherits its attributes from the C class in
Bio.Align._pairwisealigner, so they are typed by that module's stub.

Type-checked with the mypy.ini beside it, never run.
"""

from collections.abc import Callable

from typing_extensions import assert_type

from Bio.Align import PairwiseAligner

# Transitional: PairwiseAligner.__init__ itself is not annotated yet.
aligner = PairwiseAligner()  # type: ignore[no-untyped-call]

assert_type(aligner.mode, str)
aligner.mode = "local"
aligner.mode = 12345  # type: ignore[assignment]

assert_type(aligner.match_score, float | None)
aligner.match_score = 2

# The gap scores take a number or, for these three, a function.
assert_type(aligner.gap_score, float | Callable[[int, int], float])
aligner.gap_score = -1
aligner.insertion_score = lambda position, length: -2.0 * length
aligner.open_gap_score = lambda position, length: -2.0  # type: ignore[assignment]

assert_type(aligner.wildcard, str | None)
aligner.wildcard = None

assert_type(aligner.algorithm, str)
aligner.algorithm = "Gotoh"  # type: ignore[misc]
