# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Type stubs for the C extension Bio.Cluster._cluster."""

from typing import Any
from typing import overload
from typing import SupportsFloat
from typing import SupportsIndex
from typing import TypeAlias
from typing import TypeVar

from _typeshed import ReadableBuffer
from _typeshed import WriteableBuffer
from numpy.typing import NDArray

# numpy's stubs give ndarray __buffer__ only on Python 3.12+.
_Buffer: TypeAlias = ReadableBuffer | NDArray[Any]
# An array the function fills in place.
_OutBuffer: TypeAlias = WriteableBuffer | NDArray[Any]
# A distance matrix: a 1D or 2D array, or a list of rows, row i holding i
# distances.
_DistanceMatrix: TypeAlias = _Buffer | list[NDArray[Any]]

_TreeT = TypeVar("_TreeT", bound=Tree)
# Tree takes a list of any Node subclass, such as Bio.Cluster.Node.
_NodeT = TypeVar("_NodeT", bound=Node)

class Node:
    def __init__(
        self,
        left: SupportsIndex,
        right: SupportsIndex,
        distance: SupportsFloat | SupportsIndex = 0.0,
    ) -> None: ...
    @property
    def left(self) -> int: ...
    @left.setter
    def left(self, value: SupportsIndex) -> None: ...
    @property
    def right(self) -> int: ...
    @right.setter
    def right(self, value: SupportsIndex) -> None: ...
    @property
    def distance(self) -> float: ...
    @distance.setter
    def distance(self, value: SupportsFloat | SupportsIndex) -> None: ...

class Tree:
    # Only a list is accepted, not any sequence; without one, the tree is
    # empty, as Bio.Cluster.treecluster needs.
    def __new__(cls: type[_TreeT], nodes: list[_NodeT] = ..., /) -> _TreeT: ...
    def __len__(self) -> int: ...
    # Indexing gives a copy of the node, always of the base class here.
    @overload
    def __getitem__(self, index: SupportsIndex, /) -> Node: ...
    @overload
    def __getitem__(self, index: slice, /) -> list[Node]: ...
    def scale(self) -> None: ...
    # indices is an intc array of len(tree) + 1 elements, filled in place.
    def cut(self, indices: _OutBuffer, nclusters: SupportsIndex, /) -> None: ...
    def sort(self, indices: _OutBuffer, order: _Buffer, /) -> None: ...

def version() -> str: ...

# The functions below are called by the wrappers in Bio.Cluster, which
# convert each argument to a C-contiguous array of the right type first:
# double for data, weights and distances, intc for masks and cluster ids.

# Return the error and the number of times the solution was found; the
# clustering is stored in clusterid.
def kcluster(
    data: _Buffer,
    nclusters: SupportsIndex,
    mask: _Buffer,
    weight: _Buffer,
    transpose: SupportsIndex,
    npass: SupportsIndex,
    method: str,
    dist: str,
    clusterid: _OutBuffer,
    rng_seed: int | None = None,
) -> tuple[float, int]: ...
def kmedoids(
    distance: _DistanceMatrix,
    nclusters: SupportsIndex,
    npass: SupportsIndex,
    clusterid: _OutBuffer,
    rng_seed: int | None = None,
) -> tuple[float, int]: ...

# Exactly one of data and distancematrix is None; tree must be empty.
def treecluster(
    tree: Tree,
    data: _Buffer | None,
    mask: _Buffer | None,
    weight: _Buffer | None,
    transpose: SupportsIndex,
    method: str,
    dist: str,
    distancematrix: _DistanceMatrix | None,
) -> None: ...
def somcluster(
    clusterids: _OutBuffer,
    celldata: _OutBuffer,
    data: _Buffer,
    mask: _Buffer,
    weight: _Buffer,
    transpose: SupportsIndex,
    inittau: SupportsFloat | SupportsIndex,
    niter: SupportsIndex,
    dist: str,
    rng_seed: int | None = None,
) -> None: ...
def clusterdistance(
    data: _Buffer,
    mask: _Buffer,
    weight: _Buffer,
    index1: _Buffer,
    index2: _Buffer,
    method: str,
    dist: str,
    transpose: SupportsIndex,
) -> float: ...
def clustercentroids(
    data: _Buffer,
    mask: _Buffer,
    clusterid: _Buffer,
    method: str,
    transpose: SupportsIndex,
    cdata: _OutBuffer,
    cmask: _OutBuffer,
) -> None: ...

# Each row of distancematrix is filled in place.
def distancematrix(
    data: _Buffer,
    mask: _Buffer,
    weight: _Buffer,
    transpose: SupportsIndex,
    dist: str,
    distancematrix: list[NDArray[Any]],
) -> None: ...

# The C function ignores keyword arguments.
def pca(
    data: _Buffer,
    columnmean: _OutBuffer,
    coordinates: _OutBuffer,
    pc: _OutBuffer,
    eigenvalues: _OutBuffer,
    /,
) -> None: ...
