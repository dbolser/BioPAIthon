# Copyright 2000 Brad Chapman.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Extract information from alignment objects (REMOVED).

Bio.Align.AlignInfo has been removed from BioPAIthon. Its last remaining content, the SummaryInfo class, was deprecated in release 1.86; everything else was removed in 1.86. Instead of ``AlignInfo.SummaryInfo(msa).get_column(i)``, use the ``alignment`` property of the MultipleSeqAlignment object to obtain a new-style Alignment object, and take the column with ``msa.alignment[:, i]``. For consensus sequences and position-specific scores, create a ``Bio.motifs.Motif`` from that Alignment object. See DEPRECATED.rst.
"""

raise ImportError(
    "Bio.Align.AlignInfo has been removed from BioPAIthon. Its last remaining content, the SummaryInfo class, was deprecated in release 1.86; everything else was removed in 1.86. Instead of ``AlignInfo.SummaryInfo(msa).get_column(i)``, use the ``alignment`` property of the MultipleSeqAlignment object to obtain a new-style Alignment object, and take the column with ``msa.alignment[:, i]``. For consensus sequences and position-specific scores, create a ``Bio.motifs.Motif`` from that Alignment object. See DEPRECATED.rst."
)
