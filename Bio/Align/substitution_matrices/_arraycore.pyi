# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Type stubs for the C extension Bio.Align.substitution_matrices._arraycore."""

from collections.abc import Sequence
from typing import Any

from numpy.typing import NDArray

class Array(NDArray[Any]):
    def __array_finalize__(self, obj: NDArray[Any] | None, /) -> None: ...
    # None only on an array that never had an alphabet set, such as a view
    # of a plain ndarray; Bio.Align.substitution_matrices.Array always sets
    # one. The alphabet can be set only once.
    @property
    def alphabet(self) -> Sequence[Any] | Any: ...
    @alphabet.setter
    def alphabet(self, value: Sequence[Any]) -> None: ...
