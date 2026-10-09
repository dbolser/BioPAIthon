# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Type stubs for the C extension Bio.PDB.ccealign."""

from collections.abc import Sequence
from typing import Any
from typing import Final
from typing import final
from typing import SupportsFloat
from typing import SupportsIndex
from typing import TypeAlias

from _typeshed import structseq
from numpy.typing import NDArray

# Any sequence of (x, y, z) sequences; each value is read with
# PyFloat_AsDouble.
_Coordinates: TypeAlias = (
    Sequence[Sequence[SupportsFloat | SupportsIndex]] | NDArray[Any]
)

@final
class CEAlignment(structseq[Any], tuple[list[list[int]], float, int]):
    __match_args__: Final = ("path", "z_score", "length")
    # Two lists of aligned indices, into the first and the second coordinates.
    @property
    def path(self) -> list[list[int]]: ...
    @property
    def z_score(self) -> float: ...
    # The number of aligned pairs.
    @property
    def length(self) -> int: ...

def run_cealign(
    coordsA: _Coordinates,
    coordsB: _Coordinates,
    fragmentSize: SupportsIndex = 8,
    gapMax: SupportsIndex = 30,
    /,
) -> list[CEAlignment]: ...
