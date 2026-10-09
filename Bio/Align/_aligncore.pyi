# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Type stubs for the C extension Bio.Align._aligncore."""

from typing import Any
from typing import final
from typing import SupportsIndex
from typing import TypeAlias

from _typeshed import WriteableBuffer
from numpy.typing import NDArray

# numpy's stubs give ndarray __buffer__ only on Python 3.12+.
_Buffer: TypeAlias = WriteableBuffer | NDArray[Any]

@final
class PrintedAlignmentParser:
    # eol is a single byte.
    def __new__(cls, eol: bytes | bytearray = b"\n") -> PrintedAlignmentParser: ...
    # Return the number of columns read and the ungapped sequence.
    def feed(self, line: bytes, offset: SupportsIndex = 0, /) -> tuple[int, bytes]: ...
    # arr is a 2D intp array of the shape below, filled in place.
    def fill(self, arr: _Buffer, /) -> None: ...
    @property
    def shape(self) -> tuple[int, int]: ...
