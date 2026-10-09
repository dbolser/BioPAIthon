# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Type stubs for the C extension Bio.SeqIO._twoBitIO."""

from typing import Any
from typing import SupportsIndex
from typing import TypeAlias

from _typeshed import ReadableBuffer
from numpy.typing import NDArray

# numpy's stubs give ndarray __buffer__ only on Python 3.12+.
_Buffer: TypeAlias = ReadableBuffer | NDArray[Any]

# Decode the bases start:end:step from data, the packed 2bit bytes that hold
# them, applying nBlocks (as N) and maskBlocks (as lower case), each a 2D
# uint32 array of (start, end) rows. data is parsed as a read-only bytes-like
# object ("y#"), which refuses a bytearray or a memoryview.
def convert(
    data: bytes | NDArray[Any],
    start: SupportsIndex,
    end: SupportsIndex,
    step: SupportsIndex,
    nBlocks: _Buffer,
    maskBlocks: _Buffer,
) -> bytes: ...
