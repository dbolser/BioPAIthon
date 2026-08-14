# Copyright 2002 by Jeffrey Chang.
# Copyright 2016, 2019, 2020 by Markus Piotrowski.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Pairwise sequence alignment using a dynamic programming algorithm (REMOVED).

Bio.pairwise2 has been removed from BioPAIthon. It was deprecated in Biopython release 1.80 (upstream Biopython still ships it, deprecated). Please use Bio.Align.PairwiseAligner instead. Migration is not mechanical: the defaults and semantics differ. In particular, pairwise2's default gap score is 0, while PairwiseAligner's has been -1 since release 1.86, so the two can rank alignments differently unless you set match, mismatch and gap scores explicitly. See DEPRECATED.rst.
"""

raise ImportError(
    "Bio.pairwise2 has been removed from BioPAIthon. It was deprecated in Biopython release 1.80 (upstream Biopython still ships it, deprecated). Please use Bio.Align.PairwiseAligner instead. Migration is not mechanical: the defaults and semantics differ. In particular, pairwise2's default gap score is 0, while PairwiseAligner's has been -1 since release 1.86, so the two can rank alignments differently unless you set match, mismatch and gap scores explicitly. See DEPRECATED.rst."
)
