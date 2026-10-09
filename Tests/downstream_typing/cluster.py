# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""What a downstream ``mypy --strict`` run sees of Bio.Cluster.

Type-checked with the mypy.ini beside it, never run.

The arrays are not pinned to exact dtypes with assert_type, as numpy's stubs
spell those differently from release to release. Assignments to integer and
floating arrays check what kind of array each result is instead.
"""

from typing import Any

import numpy as np
from numpy.typing import NDArray
from typing_extensions import assert_type

from Bio.Cluster import kcluster
from Bio.Cluster import Node
from Bio.Cluster import pca
from Bio.Cluster import Record
from Bio.Cluster import somcluster
from Bio.Cluster import Tree
from Bio.Cluster import treecluster

data = [[1.0, 2.0], [1.5, 2.5], [8.0, 9.0]]

# k-means gives the cluster of each item, the error and how often the
# solution was found.
clusterid, error, nfound = kcluster(data, nclusters=2, rng_seed=1)
assert_type(error, float)
assert_type(nfound, int)
ids: NDArray[np.integer[Any]] = clusterid
not_floats: NDArray[np.floating[Any]] = clusterid  # type: ignore[assignment]
clusterid, error = kcluster(data)  # type: ignore[misc]

# An int works as the transpose flag, as the Tutorial uses it. Counts and
# the transpose flag can be NumPy integers.
kcluster(data, transpose=1)
kcluster(data, nclusters=np.int64(2), npass=np.int64(5), transpose=np.int64(0))
kcluster(data, method=1)  # type: ignore[arg-type]

# somcluster reads data.shape before converting data to an array, so a list
# fails at run time. The grid sizes and inittau can be NumPy scalars.
grid_ids, celldata = somcluster(
    np.array(data), nxgrid=np.int64(3), inittau=np.float32(0.05)
)
centroids: NDArray[np.floating[Any]] = celldata
somcluster(data)  # type: ignore[arg-type]

columnmean, coordinates, components, eigenvalues = pca(data)

# A Tree takes a list of Nodes, including a list of Bio.Cluster.Node.
nodes = [Node(1, 2, 0.2), Node(0, -1, 0.5)]
tree = Tree(nodes)
assert_type(tree, Tree)
assert_type(len(tree), int)
assert_type(tree[0].left, int)
assert_type(tree[-1].distance, float)
Tree(tuple(nodes))  # type: ignore[arg-type]
# A Tree has a length and indexes, but iterating one raises TypeError.
for node in tree:  # type: ignore[attr-defined]
    pass

# A slice gives copies of the nodes, to edit and build a new Tree from.
copies = tree[:]
copies[0] = Node(0, 1, 0.2)
copies[1].left = 2
copies[1].distance = 1
copies[1].right = 0.5  # type: ignore[assignment]
tree = Tree(copies)

assert_type(treecluster(data), Tree)
assert_type(treecluster(None, distancematrix=[[], [1.0]]), Tree)
cut: NDArray[np.integer[Any]] = tree.cut(np.int64(2))
order: NDArray[np.integer[Any]] = tree.sort([1.0, 2.0, 3.0])

# A Record made without a handle is filled in by hand. mask stays None
# when no value is missing, so it needs narrowing; data does not.
record = Record()
record.data = np.array(data)
record.geneid = ["a", "b", "c"]
nrows, ncolumns = record.data.shape
record.mask.shape  # type: ignore[union-attr]
record.save("jobname", record.treecluster(), record.treecluster(transpose=1))
