# Copyright 2001 Brad Chapman.
# Revisions copyright 2009-2010 by Peter Cock.
# Revisions copyright 2010 by Phillip Garland.
# All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Command line wrappers for the NCBI BLAST+ tools (REMOVED).

Bio.Blast.Applications has been removed from Biopython. It provided NcbiblastnCommandline and the other NCBI BLAST+ wrappers. It wrapped command line tools through ``Bio.Application``, which was declared obsolete in release 1.79, deprecated in 1.82, and removed in 1.86. Call the tool with the standard library ``subprocess`` module instead, for example subprocess.run(['blastn', '-query', 'in.fasta', '-db', 'nt', '-outfmt', '5'], check=True), then parse the result with Bio.Blast.parse. See DEPRECATED.rst.
"""

raise ImportError(
    "Bio.Blast.Applications has been removed from Biopython. It provided NcbiblastnCommandline and the other NCBI BLAST+ wrappers. It wrapped command line tools through ``Bio.Application``, which was declared obsolete in release 1.79, deprecated in 1.82, and removed in 1.86. Call the tool with the standard library ``subprocess`` module instead, for example subprocess.run(['blastn', '-query', 'in.fasta', '-db', 'nt', '-outfmt', '5'], check=True), then parse the result with Bio.Blast.parse. See DEPRECATED.rst."
)
