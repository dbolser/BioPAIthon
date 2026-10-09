# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Type stubs for the C extension Bio.motifs._pwm."""

from typing import Any
from typing import TypeAlias

from _typeshed import ReadableBuffer
from _typeshed import WriteableBuffer
from numpy.typing import NDArray

# numpy's stubs give ndarray __buffer__ only on Python 3.12+.
_Buffer: TypeAlias = ReadableBuffer | NDArray[Any]
_WriteableBuffer: TypeAlias = WriteableBuffer | NDArray[Any]

# Fill scores, a 1D float32 array of len(sequence) - len(matrix) + 1 items,
# with the score of matrix at each position of sequence; NaN where a letter is
# none of A, C, G and T. matrix is a 2D float64 array with a column for each of
# those letters. sequence is parsed as a read-only bytes-like object ("y#"),
# which refuses a bytearray or a memoryview.
def calculate(
    sequence: bytes | NDArray[Any], matrix: _Buffer, scores: _WriteableBuffer
) -> None: ...
