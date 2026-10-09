# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Type stubs for the C extension Bio.PDB._bcif_helper."""

from typing import Any
from typing import TypeAlias

from _typeshed import ReadableBuffer
from _typeshed import WriteableBuffer
from numpy.typing import NDArray

# numpy's stubs give ndarray __buffer__ only on Python 3.12+.
_Buffer: TypeAlias = ReadableBuffer | NDArray[Any]
_WriteableBuffer: TypeAlias = WriteableBuffer | NDArray[Any]

# Undo BinaryCIF's IntegerPacking. packed is a 1D array of int8 or int16,
# filling out, a 1D int32 array; or of uint8 or uint16, filling a uint32 one.
# Both are in native byte order, and out must be exactly as long as the
# unpacked data.
def integer_unpack(packed: _Buffer, out: _WriteableBuffer, /) -> None: ...
