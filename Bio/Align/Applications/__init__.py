# Copyright 2009 by Peter Cock & Cymon J. Cox.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Command line wrappers for multiple sequence alignment tools (REMOVED).

Bio.Align.Applications has been removed from Biopython. It provided MuscleCommandline, MafftCommandline, ClustalOmegaCommandline and the other alignment wrappers. It wrapped command line tools through ``Bio.Application``, which was declared obsolete in release 1.79, deprecated in 1.82, and removed in 1.86. Call the tool with the standard library ``subprocess`` module instead, then read the result with Bio.Align.read or Bio.AlignIO.read. See DEPRECATED.rst.
"""

raise ImportError(
    "Bio.Align.Applications has been removed from Biopython. It provided MuscleCommandline, MafftCommandline, ClustalOmegaCommandline and the other alignment wrappers. It wrapped command line tools through ``Bio.Application``, which was declared obsolete in release 1.79, deprecated in 1.82, and removed in 1.86. Call the tool with the standard library ``subprocess`` module instead, then read the result with Bio.Align.read or Bio.AlignIO.read. See DEPRECATED.rst."
)
