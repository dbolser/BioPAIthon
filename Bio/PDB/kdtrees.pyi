# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Type stubs for the C extension Bio.PDB.kdtrees."""

from typing import Any
from typing import final
from typing import SupportsIndex
from typing import TypeAlias

from _typeshed import ReadableBuffer
from numpy.typing import NDArray
from typing_extensions import Self

# numpy's stubs give ndarray __buffer__ only on Python 3.12+.
_Buffer: TypeAlias = ReadableBuffer | NDArray[Any]

class KDTree:
    # coords is an Nx3 float64 array. Both arguments are positional only: the
    # C code silently ignores keyword arguments, so bucket_size=10 gives 1.
    def __new__(cls, coords: _Buffer, bucket_size: SupportsIndex = 1, /) -> Self: ...
    # center is a float64 array of three coordinates.
    def search(self, center: _Buffer, radius: float, /) -> list[Point]: ...
    def neighbor_search(self, radius: float, /) -> list[Neighbor]: ...
    def neighbor_simple_search(self, radius: float, /) -> list[Neighbor]: ...

@final
class Point:
    def __init__(self, index: SupportsIndex, radius: float = 0.0) -> None: ...
    @property
    def index(self) -> int: ...
    @property
    def radius(self) -> float: ...

@final
class Neighbor:
    def __init__(
        self, index1: SupportsIndex, index2: SupportsIndex, radius: float = 0.0
    ) -> None: ...
    @property
    def index1(self) -> int: ...
    @property
    def index2(self) -> int: ...
    @property
    def radius(self) -> float: ...
