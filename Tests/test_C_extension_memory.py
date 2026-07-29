# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Check that the C extensions do not grow their memory use per call.

See memory_growth.py for why this is measured by repetition rather than by
LeakSanitizer. Bio.PDB.ccealign is covered by test_PDB_CEAligner instead,
alongside the reference-ownership test that belongs with it.
"""

import unittest

from memory_growth import assert_bounded_growth
from memory_growth import requires_growth_measurement

try:
    import numpy as np
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use the C extensions."
    ) from None


@requires_growth_measurement
class AlignmentExtensionTests(unittest.TestCase):
    """Bio.Align._aligncore and Bio.Align._pairwisealigner."""

    def test_printed_alignment_parser(self):
        """Parse a printed alignment repeatedly."""
        from Bio.Align import Alignment

        lines = [b"ACGT--ACGTACGT--AC", b"ACGTTTACGT--ACGTAC"]

        def parse():
            Alignment.parse_printed_alignment(lines)

        assert_bounded_growth(self, parse)

    def test_pairwise_aligner_score(self):
        """Score a pair of sequences repeatedly."""
        from Bio.Align import PairwiseAligner

        aligner = PairwiseAligner()
        target = "AGAACTTATCGCTTGACGTAAC" * 2
        query = "AGAACTATCGCTGACGTTAAC" * 2

        def score():
            aligner.score(target, query)

        assert_bounded_growth(self, score)

    def test_pairwise_aligner_alignments(self):
        """Build and read back alignments repeatedly.

        This is the path that allocates most: the aligner holds the traceback
        matrix until the alignments iterator is exhausted.
        """
        from Bio.Align import PairwiseAligner

        aligner = PairwiseAligner()
        target = "AGAACTTATCGCTTGACGTAAC"
        query = "AGAACTATCGCTGACGTTAAC"

        def align():
            alignments = aligner.align(target, query)
            for index, alignment in enumerate(alignments):
                if index == 2:
                    break
                str(alignment)

        assert_bounded_growth(self, align)


@requires_growth_measurement
class StructureExtensionTests(unittest.TestCase):
    """Bio.PDB.kdtrees."""

    def test_kdtree_search(self):
        """Build a tree and search it repeatedly."""
        from Bio.PDB import kdtrees

        # The extension requires doubles, and coordinates within +/- 1e6.
        rng = np.random.default_rng(0)
        coords = rng.random((500, 3)) * 20.0
        center = np.array([10.0, 10.0, 10.0])

        def search():
            tree = kdtrees.KDTree(coords, 10)
            tree.search(center, 5.0)
            tree.neighbor_search(3.0)

        assert_bounded_growth(self, search)


@requires_growth_measurement
class ClusterExtensionTests(unittest.TestCase):
    """Bio.Cluster._cluster."""

    def test_kcluster_and_distancematrix(self):
        """Cluster a small matrix repeatedly."""
        from Bio import Cluster

        rng = np.random.default_rng(0)
        data = rng.random((30, 8))

        def cluster():
            Cluster.kcluster(data, nclusters=3, npass=1)
            Cluster.distancematrix(data)

        assert_bounded_growth(self, cluster)

    def test_treecluster(self):
        """Build a hierarchical clustering tree repeatedly."""
        from Bio import Cluster

        rng = np.random.default_rng(0)
        data = rng.random((30, 8))

        def cluster():
            Cluster.treecluster(data)

        assert_bounded_growth(self, cluster)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
