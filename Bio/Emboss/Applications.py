# Copyright 2001-2009 Brad Chapman.
# Revisions copyright 2009-2016 by Peter Cock.
# Revisions copyright 2009 by David Winter.
# Revisions copyright 2009-2010 by Leighton Pritchard.
# All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Command line wrappers for the EMBOSS suite (REMOVED).

Bio.Emboss.Applications has been removed from Biopython. It provided WaterCommandline, NeedleCommandline, SeqretCommandline and the other EMBOSS wrappers. It wrapped command line tools through ``Bio.Application``, which was declared obsolete in release 1.79, deprecated in 1.82, and removed in 1.86. Call the tool with the standard library ``subprocess`` module instead. Bio.Emboss.Primer3 and Bio.Emboss.PrimerSearch, which parse EMBOSS output rather than run it, are unaffected and remain available. See DEPRECATED.rst.
"""

raise ImportError(
    "Bio.Emboss.Applications has been removed from Biopython. It provided WaterCommandline, NeedleCommandline, SeqretCommandline and the other EMBOSS wrappers. It wrapped command line tools through ``Bio.Application``, which was declared obsolete in release 1.79, deprecated in 1.82, and removed in 1.86. Call the tool with the standard library ``subprocess`` module instead. Bio.Emboss.Primer3 and Bio.Emboss.PrimerSearch, which parse EMBOSS output rather than run it, are unaffected and remain available. See DEPRECATED.rst."
)
